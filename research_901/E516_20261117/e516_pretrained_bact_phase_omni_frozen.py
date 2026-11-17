#!/usr/bin/env python3
"""901/E516: no-tuning third-party model on original MoMA-*independent* phase TIFFs.

Strict evidence boundary: published pretrained Omnipose 1.1.4 bact_phase_omni
inference on source phase pixels is independent of MoMA candidate masks, but
NOT independent experimental ground truth. Controls are software-agreement
and ordinary windows, NOT labeled biological positive/negative truth.

Model not installed/weights unavailable: fail loudly, save reason separately.
E510 source pixel assets and all pre-existing E508-E515 files remain immutable.
"""
import argparse, collections, csv, hashlib, importlib.metadata, io, json, os, pathlib, platform, sys, traceback, zipfile
import numpy as np
from scipy.ndimage import find_objects

P=argparse.ArgumentParser()
P.add_argument('--gl',choices=['Pos5_GL35','Pos12_GL22'],required=True)
P.add_argument('--e510-zip',required=True)
P.add_argument('--output',required=True)
P.add_argument('--inspect-source-only',action='store_true')
P.add_argument('--preflight-count',type=int,default=8,
               help='First N unique frames to segment in minimal setup test; 0=all')
a=P.parse_args()
OUT=pathlib.Path(a.output)/a.gl;OUT.mkdir(parents=True,exist_ok=True)
FROZEN_SOURCE={
 'Pos5_GL35':'745815a73bc1e2eaa01a2cf99f4edd9c530064d10c22bc5421f6ff5189629d8d',
 'Pos12_GL22':'17f80af916f376dfa380b30e8a8b2e8459358a4b66e13e21e23f637aaa08cfa9'
}
MODEL_ID='bact_phase_omni'; MODEL_RELEASE='omnipose==1.1.4'; NO_FINETUNE=True
INFERENCE={'gpu':False,'model_type':MODEL_ID,'net_avg':False,
    'channels':[0,0],'diameter':None,'normalize':True,
    'invert':False,'omni':True,'flow_threshold':0.0,
    'mask_threshold':0.0,'min_size':15,'tile':False,
    'resample':True,'augment':False}
log={'stage':'E516_3RD_PARTY_MODEL_NO_TUNING',
     'model_requested':MODEL_ID,'pypi_stable_version':MODEL_RELEASE,
     'no_training':True,'no_MoMA_source_candidates_or_labels_in_model_input':True,
     'not_external_human_segmentation_ground_truth':True,
     'original_pixel_zero_based_frame':'TIFF_time=MoMA_frame; E510 labels were TIFF_time+1',
     'source_gl':a.gl,'software_controls_not_biological_truth':True,
     'parameter_frozen_before_model_inference':INFERENCE,'status':'PREPARED'}
def write_log():
 (OUT/'E516_pretrained_model_source_provenance_and_status.json').write_text(
     json.dumps(log,indent=2,ensure_ascii=False),encoding='utf8')
write_log()
try:
    in_bytes=pathlib.Path(a.e510_zip).read_bytes()
    digest=hashlib.sha256(in_bytes).hexdigest()
    assert digest==FROZEN_SOURCE[a.gl],('E510_SOURCE_SHA_MISMATCH',digest)
    with zipfile.ZipFile(io.BytesIO(in_bytes)) as z:
        assert z.testzip() is None
        names=z.namelist()
        source_key=[n for n in names if n.endswith('E510_sealed_event_key_do_not_open_for_blind_review.json')]
        assert len(source_key)==1
        frame_arrays={}
        provenance=collections.defaultdict(list)
        for name in names:
            if not name.endswith('_exact_float32_phase_frames.npz'):continue
            with np.load(io.BytesIO(z.read(name)),allow_pickle=False) as x:
                image=x['phase_pixels'];labels=x['frames'];old_offsets=x['offsets']
                assert image.dtype==np.float32 and image.shape==(9,300,106)
                assert old_offsets.tolist()==[-3,-2,-1,0,1,2,3,5,8]
                for lab,im in zip(labels,image):
                    t=int(lab)-1
                    if t in frame_arrays:
                        assert np.array_equal(frame_arrays[t],im),(a.gl,t,'DUPLICATED_PIXEL_MISMATCH')
                    else:
                        frame_arrays[t]=np.ascontiguousarray(im)
                    provenance[t].append(name)
        expected=119 if a.gl=='Pos5_GL35' else 84
        assert len(frame_arrays)==expected,(len(frame_arrays),expected)
        sha_by_t={str(t):hashlib.sha256(a.tobytes()).hexdigest() for t,a in frame_arrays.items()}
        (OUT/'E516_frozen_original_pixel_sha256.json').write_text(json.dumps({
            'gl':a.gl,'TIFF_source_zero_based_TIFF_frames':sha_by_t,
            'unblinded_software_event_class_not_used':True,
            'source_attachment_E510_zip_SHA256':digest},indent=2),encoding='utf8')
        log.update({'source_e510_artifact_sha256':digest,
                    'source_zip_member_crc32':'ALL_PASS',
                    'unique_source_TIFF_frames':len(frame_arrays),
                    'raw_phase_image_shape':[300,106],
                    'raw_phase_dtype':'float32',
                    'source_data_authenticity':'PASS'})
        write_log()
        if a.inspect_source_only:
            log.update({'status':'SOURCE_INTEGRITY_PASS_MODEL_NOT_RUN'})
            write_log()
            print('E516_SOURCE_INPUT_PASS',a.gl,len(frame_arrays),flush=True)
            sys.exit(0)
        # The sealed E510 key is not read until after inference and masks have
        # been saved and hashed; named windows carry no class at inference.
        sealed_key=z.read(source_key[0])
    # Actual model must be the pinned, non-yanked published stable release.
    installed=importlib.metadata.version('omnipose')
    assert installed=='1.1.4',('MODEL_VERSION_NOT_FROZEN',installed)
    from cellpose_omni import models
    import torch
    log.update({'installed_omnipose_version':installed,
                'python_version':sys.version,'torch_version':torch.__version__,
                'torch_cpu_only':True,
                'platform':platform.platform()})
    # 1.1.4 automatically maps bact_phase_omni to nchan=2 and nclasses=3 abstract;
    # actual torch head output is 4 for 2D.
    model=models.CellposeModel(gpu=False,model_type=MODEL_ID,net_avg=False,
        nclasses=3,nchan=2,dim=2,omni=True)
    weight_files=[pathlib.Path(p) for p in model.pretrained_model]
    assert len(weight_files)==1 and weight_files[0].is_file()
    weight=weight_files[0]
    assert weight.stat().st_size>1_000_000
    model_sha=hashlib.sha256(weight.read_bytes()).hexdigest()
    log.update({'weight_filename':weight.name,'weight_bytes':weight.stat().st_size,
                'weight_sha256':model_sha,'real_pretrained_model_weight_loaded':True,
                'nchan':model.nchan,'nclasses_output':model.nclasses})
    write_log()
    all_times=sorted(frame_arrays)
    test_times=all_times[:a.preflight_count] if a.preflight_count>0 else all_times
    # All model calls are made on exactly the original single-channel float TIFF,
    # never on 8-bit PNG rescalings, MoMA segmentation masks, probabilities,
    # software candidate locations, or event labels.
    masks={};features=[]
    for t in test_times:
        x=frame_arrays[t]
        out=model.eval(x,channels=INFERENCE['channels'],diameter=INFERENCE['diameter'],
            normalize=INFERENCE['normalize'],invert=INFERENCE['invert'],
            omni=INFERENCE['omni'],flow_threshold=INFERENCE['flow_threshold'],
            mask_threshold=INFERENCE['mask_threshold'],min_size=INFERENCE['min_size'],
            tile=INFERENCE['tile'],resample=INFERENCE['resample'],
            augment=INFERENCE['augment'],net_avg=False,show_progress=False)
        assert isinstance(out,tuple) and len(out)>=1
        m=np.asarray(out[0]);m=np.squeeze(m)
        assert m.shape==x.shape and np.issubdtype(m.dtype,np.integer),(m.shape,m.dtype)
        assert m.min()>=0
        m=m.astype(np.uint16)
        masks[t]=m
        ids,counts=np.unique(m,return_counts=True)
        nonzero=counts[ids>0]
        features.append({'t_source_original_TIFF_zero_based':t,
            'n_instances':len(nonzero),'foreground_fraction':float((m>0).mean()),
            'median_instance_area_px':float(np.median(nonzero)) if len(nonzero) else None,
            'source_phase_sha256':sha_by_t[str(t)],
            'inference_mask_sha256':hashlib.sha256(m.tobytes()).hexdigest()})
        print('E516_TRUE_PRETRAINED_MASK',a.gl,t,len(nonzero),flush=True)
    np.savez_compressed(OUT/'E516_original_TIFF_pretrained_segmentation_masks.npz',
                        time_original_TIFF_zero_based=np.array(test_times,dtype='int32'),
                        labels=np.stack([masks[t] for t in test_times]),
                        source_pixel_sha256=np.array([sha_by_t[str(t)] for t in test_times],dtype='U64'))
    (OUT/'E516_pretrained_no_tuning_frame_qc.json').write_text(
      json.dumps({'model':MODEL_ID,'weights_SHA256':model_sha,'frame_QC':features},indent=2),encoding='utf8')
    digest_out=hashlib.sha256((OUT/'E516_original_TIFF_pretrained_segmentation_masks.npz').read_bytes()).hexdigest()
    log.update({'status':'PRETRAINED_INFERENCE_COMPLETE_NO_BIO_TRUTH',
        'model_inferred_original_unique_frames':len(test_times),
        'mask_archive_sha256':digest_out,
        'min_mask_count':min(x['n_instances'] for x in features),
        'median_mask_count':float(np.median([x['n_instances'] for x in features])),
        'no_event_classification_or_accuracy_claimed':True})
    write_log()
    print('E516_VERIFIED_PRETRAINED_INFERENCE_PASS',a.gl,MODEL_ID,len(test_times),model_sha,flush=True)
except SystemExit:
    raise
except BaseException as e:
    log.update({'status':'FAILED_NO_INDEPENDENT_MODEL_INFERENCE',
                'exception_type':type(e).__name__,'exception_message':str(e),
                'traceback_tail':traceback.format_exc()[-4500:]})
    write_log()
    raise