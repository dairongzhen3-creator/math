#!/usr/bin/env python3
"""901 E554 Colab recovery without Drive filesystem mount.

Only operates on completed source-certified /content output. Makes a user-downloadable
backup first, then (if Colab OAuth is available) uploads to the two verified Drive
roots with per-file MD5 and full-byte SHA256 API readback. Never overwrites originals.
"""
import datetime
import hashlib
import io
import json
from pathlib import Path
import zipfile

OUT = Path('/content/901_E554_official_extraction')
PUBLIC_PARENT = '1KKGKGvDkpTxVpjfg5lD38mKOW8er5npb'  # 03_工具与公共数据
PROJECT_PARENT = '1vs3ac-hYrKlt0zYCcbZgQGfLaAjBvwYI' # 02_项目大数据/901_过程数据
PUBLIC_FOLDER = 'S-BIAD2506_MoMA_E554_original_small_sidecars'
PROJECT_FOLDER = 'E554_两关键帧_源记录_20261118'


def file_hash(path, algo):
    digest = hashlib.new(algo)
    with Path(path).open('rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def source_is_public(name):
    return (name.startswith(('Pos5_GL35__', 'Pos12_GL22__')) and (
        name.endswith(('moma.log', 'mm.properties', 'curation.moma', '.csv'))
        or 'curation.moma__BKP' in name))


def make_backup(files, provenance):
    dest = Path('/content/901_E554_official_source_backup.zip')
    manifest = {'study': '901', 'stage': 'E554', 'certification': 'OFFICIAL_SOURCE_PASS',
                'files': provenance, 'source_result': 'E554_full_source_manifest.json'}
    with zipfile.ZipFile(dest, mode='w', compression=zipfile.ZIP_DEFLATED, compresslevel=6,
                         allowZip64=True) as archive:
        archive.writestr('E554_rescue_files_sha256.json',
                         json.dumps(manifest, ensure_ascii=False, sort_keys=True, indent=2)+'\n')
        for src in files:
            root = '03_public_originals' if source_is_public(src.name) else '02_project_frames_and_proofs'
            archive.write(src, arcname=f'{root}/{src.name}')
    # The file contents and manifest have been generated in the current Colab cloud runtime.
    print('BACKUP_ZIP_CREATED', str(dest), 'bytes', dest.stat().st_size, 'SHA256', file_hash(dest, 'sha256'))
    return dest


def ensure_folder(svc, parent, name):
    q = f"'{parent}' in parents and name='{name.replace(chr(39), chr(92)+chr(39))}' and trashed = false"
    items = svc.files().list(q=q, spaces='drive', fields='nextPageToken,files(id,name,mimeType)',
                             pageSize=100).execute().get('files', [])
    assert len(items) <= 1, ('AMBIGUOUS_DUPLICATE_FOLDER', name, items)
    if items:
        assert items[0]['mimeType'] == 'application/vnd.google-apps.folder'
        return items[0]['id']
    result = svc.files().create(body={'name': name, 'parents': [parent],
                                     'mimeType': 'application/vnd.google-apps.folder'},
                                fields='id,name,mimeType').execute()
    return result['id']


def remote_sha256(svc, fid):
    from googleapiclient.http import MediaIoBaseDownload
    buf = io.BytesIO()
    req = svc.files().get_media(fileId=fid)
    downloader = MediaIoBaseDownload(buf, req, chunksize=4 << 20)
    done = False
    while not done:
        _, done = downloader.next_chunk(num_retries=3)
    return hashlib.sha256(buf.getvalue()).hexdigest()


def upload_no_replace(svc, parent, path, known_sha):
    from googleapiclient.http import MediaFileUpload
    safe_name = path.name.replace("'", "\\'")
    q = f"'{parent}' in parents and name='{safe_name}' and trashed = false"
    items = svc.files().list(q=q, spaces='drive', fields='nextPageToken,files(id,name,size,md5Checksum)',
                             pageSize=100).execute().get('files', [])
    assert len(items) <= 1, ('AMBIGUOUS_DUPLICATE', parent, path.name)
    if items:
        doc = items[0]
        assert int(doc.get('size', -1)) == path.stat().st_size, ('EXISTING_SIZE_CONFLICT',path.name)
        assert doc.get('md5Checksum') == file_hash(path, 'md5'), ('EXISTING_MD5_CONFLICT', path.name)
        fid = doc['id']
        status = 'PRESERVED_IDENTICAL'
    else:
        media = MediaFileUpload(str(path), mimetype='application/octet-stream', resumable=True,
                                chunksize=4 << 20)
        doc = svc.files().create(body={'name':path.name,'parents':[parent]},
                                 media_body=media, fields='id,name,size,md5Checksum').execute()
        assert int(doc.get('size', -1)) == path.stat().st_size, ('UPLOAD_SIZE_MISMATCH',path.name)
        assert doc.get('md5Checksum') == file_hash(path, 'md5'), ('UPLOAD_MD5_MISMATCH',path.name)
        fid = doc['id']
        status = 'CREATED'
    remote_sha = remote_sha256(svc, fid)
    assert remote_sha == known_sha, ('INDEPENDENT_DRIVE_API_SHA256_MISMATCH', path.name)
    return {'name':path.name, 'file_id':fid,'size':path.stat().st_size,
            'sha256':known_sha,'readback_sha256':remote_sha,'status':status+'_REMOTE_SHA_PASS'}


def rescue():
    if not OUT.is_dir() or not (OUT/'E554_full_source_manifest.json').is_file():
        print('NO_PERSISTENT_SOURCE_OUTPUT_IN_THIS_COLAB_RUNTIME')
        print('Please rerun the first two code cells of the original E554 notebook, then run this rescue cell.')
        return False
    manifest = json.loads((OUT/'E554_full_source_manifest.json').read_text(encoding='utf-8'))
    assert manifest.get('status') == 'OFFICIAL_ZIP64_CENTRAL_SHA_CRC_CURATIONS_AND_SOURCE_IMAGE_FRAMES_PASS', manifest.get('status')
    assert {g.get('GL') for g in manifest.get('groups',[])} == {'Pos5_GL35','Pos12_GL22'}
    files = sorted(p for p in OUT.iterdir() if p.is_file())
    assert len(files) >= 60, ('MISSING_SOURCE_OUTPUT_FILES', len(files))
    assert all(p.name == p.name.split('/')[-1] for p in files)
    provenance = [{'name':p.name, 'bytes':p.stat().st_size, 'sha256':file_hash(p,'sha256'),
                   'target':'03_public' if source_is_public(p.name) else '02_project'} for p in files]
    print('OFFICIAL_SOURCE_OUTPUT_PRESENT',len(files),'files')
    backup = make_backup(files, provenance) # Always make backup before attempting auth
    try:
        import google.colab.auth
        google.colab.auth.authenticate_user()
        import google.auth
        from googleapiclient.discovery import build
        creds,_ = google.auth.default(scopes=['https://www.googleapis.com/auth/drive'])
        svc = build('drive','v3',credentials=creds, cache_discovery=False)
        for id_, expected in [(PUBLIC_PARENT,'03_工具与公共数据'),(PROJECT_PARENT,'901_过程数据')]:
            parent = svc.files().get(fileId=id_,fields='id,name,mimeType,trashed').execute()
            assert parent['name'] == expected and parent['mimeType']=='application/vnd.google-apps.folder' and not parent.get('trashed'), ('WRONG_DRIVE_DESTINATION',parent)
        public_folder = ensure_folder(svc, PUBLIC_PARENT, PUBLIC_FOLDER)
        project_folder = ensure_folder(svc, PROJECT_PARENT, PROJECT_FOLDER)
        receipts = []
        for i,record in enumerate(provenance,1):
            src=OUT/record['name']
            folder=public_folder if record['target']=='03_public' else project_folder
            receipts.append(upload_no_replace(svc,folder,src,record['sha256']))
            if i % 8 == 0 or i == len(provenance):
                print('DRIVE_API_SHA_ROUNDTRIP',i,'/',len(provenance),flush=True)
        result = {'study':'901','stage':'E554','source_status':manifest['status'],
                  'Drive_API_readback':'ALL_FILES_SHA256_VERIFIED',
                  'public_folder_id':public_folder,'project_folder_id':project_folder,
                  'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  'biological_division_truth':'UNRESOLVED','files':receipts}
        receipt=OUT/'E554_Colab_DriveAPI_SHA_readback.json'
        receipt.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf-8')
        upload_no_replace(svc,project_folder,receipt,file_hash(receipt,'sha256'))
        print('SOURCE_INTAKE_SUCCESS_AND_DRIVE_API_SHA_ROUNDTRIP_PASS',len(receipts),'files')
        print('Drive 03 public ID', public_folder)
        print('Drive 02 901 project ID', project_folder)
        return True
    except Exception as exc:
        print('DRIVE_API_AUTH_OR_TRANSFER_FAILED',type(exc).__name__,str(exc)[:1200])
        print('EXTRACTION_UNCHANGED. BACKUP_ARCHIVE_READY', backup)
        print('Your extraction remains in Colab while this runtime is alive.')
        print('Browser download will be offered; you can upload that ZIP to this ChatGPT conversation for cloud archival.')
        try:
            from google.colab import files as colab_files
            colab_files.download(str(backup))
        except Exception as download_error:
            print('Browser download could not be initiated:',type(download_error).__name__,str(download_error)[:500])
        return False

if __name__=='__main__':
    rescue()