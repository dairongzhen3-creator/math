#!/usr/bin/env python3
"""Index public catalogue metadata; parsing never counts as reading or proof review.

Local: python our_review/build_review_index.py --root .
GitHub runner: --remote --publish (GITHUB_TOKEN and GITHUB_REPOSITORY required).
Only generated paths under our_review/ are published. No manuscript code is run.
"""
from __future__ import annotations
import argparse
import base64
from collections import Counter
import csv
import hashlib
import html
import io
import json
import os
from pathlib import Path
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

BASELINE = 'adc7f1241b42e322a6451854ab7e4b4c146bf78a'
CATALOGUE_BLOB = '8148d2c9077e611d5665a2ded93da68e5b608f0f'
ALLOWED_REPOSITORY = 'dairongzhen3-creator/math'
STATUSES = {'UNREAD', 'ABSTRACT_READ', 'PARTIAL_TEXT_READ', 'FULL_TEXT_READ', 'PROOF_CHECKED_PARTIAL'}
TAGS = {
    'quantum-information': ['quantum channel', 'amplitude damping', 'dephasing', 'lindblad', 'entanglement', 'quantum capacity', 'von neumann entropy'],
    'spectral-inverse': ['spectral gap', 'inverse problem', 'reconstruction', 'hankel', 'prony', 'autocorrelation', 'relaxation', 'spectral measure'],
    'critical-phenomena': ['ising', 'potts', 'spin glass', 'critical', 'phase transition', 'xy model', 'heisenberg'],
    'information-learning': ['information', 'learning', 'memory', 'sample complexity', 'compression', 'regression', 'communication'],
    'modular-conformal': ['theta', 'modular form', 'vertex operator', 'conformal', 'verlinde', 'ramanujan', 'nonunitary'],
    'gauge-field': ['gauge', 'yang-mills', 'gauss law', "gauss's law", 'particle', 'mass gap', 'coulomb', 'maxwell'],
    'stochastic-dynamics': ['mixing', 'ergodic', 'markov', 'hydrodynamic', 'stochastic', 'quench', 'dynamical'],
    'geometry-topology': ['topology', 'homology', 'homotopy', 'manifold', 'curvature', 'betti', 'metric geometry'],
    'computation-certificates': ['algorithm', 'numerical', 'semidefinite', 'polynomial-time', 'approximation', 'certificate'],
    'pde-transport': ['navier', 'vlasov', 'schrodinger', 'schrödinger', 'reaction-diffusion', 'transport', 'dissipation', 'kinetic']
}

def git_blob(data: bytes) -> str:
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()

def clean(text: str) -> str:
    text = re.sub(r'<[^>]+>', '', text)
    return ' '.join(html.unescape(text).split())

def table(text: str) -> str:
    return clean(text).replace('|', '&#124;')

def parse_catalogue(text: str) -> list[dict]:
    families = list(re.finditer(r'^\*\*(\d{3})\.\s*(.*?)\*\*', text, re.M))
    rows = []
    for i, match in enumerate(families):
        stop = families[i + 1].start() if i + 1 < len(families) else len(text)
        block = text[match.end():stop]
        papers = list(re.finditer(r'&emsp;\[(.*?)\]\((preprints/[^\n]+?\.pdf)\)', block))
        description = clean(block[:papers[0].start()]) if papers else ''
        for j, paper in enumerate(papers):
            end = papers[j + 1].start() if j + 1 < len(papers) else len(block)
            abstract = clean(block[paper.end():end])
            path = paper.group(2)
            title = clean(paper.group(1))
            searchable = (match.group(2) + ' ' + title + ' ' + abstract).lower()
            tags = [tag for tag, terms in TAGS.items() if any(term in searchable for term in terms)]
            rows.append({
                'paper_id': f'F{match.group(1)}-P{j+1:02d}',
                'family_id': match.group(1),
                'family_title': clean(match.group(2)),
                'family_description': description,
                'title': title, 'paper_path': path,
                'abstract': abstract,
                'contents_source_line': text.count('\n', 0, match.end() + paper.start()) + 1,
                'candidate_tags_automatic': tags,
                'tagging_is_semantic_review': False,
                'formalization_link_listed_for_family': bool(re.search(r'lean/docs/\d+\.md', description)),
            })
    return rows

def parse_disciplines(text: str) -> dict[str, str]:
    discipline = 'UNASSIGNED'
    result = {}
    for line in text.splitlines():
        section = re.match(r'\\cataloguesection\{([^{}]+)\}', line)
        if section:
            discipline = section.group(1)
        entry = re.match(r'\\resultentry\{(\d{3})\}', line)
        if entry:
            result[entry.group(1)] = discipline
    return result

def render(catalogue: bytes, overview: bytes, state: dict) -> dict[str, str]:
    if git_blob(catalogue) != CATALOGUE_BLOB:
        raise ValueError('Catalogue blob differs from fixed baseline; do not silently change census.')
    rows = parse_catalogue(catalogue.decode('utf-8'))
    families = {row['family_id'] for row in rows}
    paths = [row['paper_path'] for row in rows]
    if len(rows) != 722 or len(families) != 372 or len(set(paths)) != 722:
        raise ValueError(f'Census mismatch: entries={len(rows)}, families={len(families)}, unique_paths={len(set(paths))}')
    if families != {f'{i:03d}' for i in range(1, 373)}:
        raise ValueError('Family identifiers are not the complete 001-372 sequence.')
    if any(not row['abstract'] for row in rows):
        raise ValueError('Empty abstract encountered.')
    if state.get('baseline_commit') != BASELINE:
        raise ValueError('Reading state has a different baseline.')
    paper_state = state.get('papers', {})
    unknown = set(paper_state) - set(paths)
    if unknown:
        raise ValueError(f'Reading state contains unknown paths: {sorted(unknown)}')
    read_families = set(state.get('abstract_read_families', []))
    if not read_families <= families:
        raise ValueError('Unknown abstract-read family.')
    disciplines = parse_disciplines(overview.decode('utf-8'))
    if set(disciplines) != families:
        raise ValueError('Overview discipline map does not cover exactly the catalogue families.')
    for row in rows:
        entry = paper_state.get(row['paper_path'], {})
        default = 'ABSTRACT_READ' if row['family_id'] in read_families else 'UNREAD'
        status = entry.get('status', default)
        if status not in STATUSES:
            raise ValueError(f'Invalid reading status: {status}')
        row.update({'discipline': disciplines[row['family_id']], 'reading_status': status,
                    'read_scope': entry.get('read_scope', []), 'review_date': entry.get('date'),
                    'proof_verified': bool(entry.get('proof_verified', False)),
                    'review_note': entry.get('note', ''),
                    'baseline_commit': BASELINE,
                    'source_url': f'https://github.com/openai/math/blob/{BASELINE}/' + row['paper_path']})
    counts = Counter(row['reading_status'] for row in rows)
    audit = {'baseline_commit': BASELINE, 'catalogue_git_blob': git_blob(catalogue),
             'catalogue_bytes': len(catalogue), 'catalogue_sha256': hashlib.sha256(catalogue).hexdigest(),
             'overview_git_blob': git_blob(overview),
             'manuscripts': len(rows), 'unique_pdf_paths': len(set(paths)), 'families': len(families),
             'nonempty_abstracts': sum(bool(row['abstract']) for row in rows),
             'reading_status_counts': dict(sorted(counts.items())),
             'proof_verified_count': sum(row['proof_verified'] for row in rows),
             'discipline_manuscript_counts': dict(sorted(Counter(row['discipline'] for row in rows).items())),
             'automatic_tags_are_not_reading': True,
             'individual_pdf_bytes_checked': False,
             'original_manuscripts_modified': False}
    output = {'our_review/index_audit.json': json.dumps(audit, ensure_ascii=False, indent=2) + '\n'}
    output['our_review/catalogue_rows.jsonl'] = ''.join(json.dumps(row, ensure_ascii=False) + '\n' for row in rows)
    header = '# 722 篇逐篇审读台账\n\n'
    header += f'固定基线 `{BASELINE}`；脚本逐项解析 {len(rows)} 篇、{len(families)} 家族、{len(set(paths))} 个不同 PDF 路径。\n\n'
    header += '**登记不等于阅读。** 状态仅来自 review_state.json；自动标签不表示语义审查或无关判决。原文结果均为作者陈述，未由此脚本验证。\n\n'
    header += '状态计数：' + '; '.join(f'{key}={value}' for key, value in sorted(counts.items())) + '\n\n'
    header += '| ID | 手稿（原文链接） | 学科 | 自动候选标签 | 实际阅读状态 |\n|---|---|---|---|---|\n'
    lines = [f"| {r['paper_id']} | [{table(r['title'])}]({r['source_url']}) | {table(r['discipline'])} | {', '.join(r['candidate_tags_automatic']) or 'NEEDS_MANUAL_SCREEN'} | {r['reading_status']} |" for r in rows]
    output['our_review/00_722_master_index.md'] = header + '\n'.join(lines) + '\n'
    csv_buffer = io.StringIO(newline='')
    fields = ['paper_id', 'family_id', 'family_title', 'title', 'discipline', 'paper_path', 'reading_status', 'proof_verified', 'contents_source_line', 'source_url']
    writer = csv.DictWriter(csv_buffer, fieldnames=fields, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows)
    output['our_review/722_master_index.csv'] = '\ufeff' + csv_buffer.getvalue()
    family_lines = ['# 372 家族摘要阅读包\n\n各文件只含官方家族摘要与论文列表。取得文件不计作已读。\n']
    family_rows = {fid: [r for r in rows if r['family_id'] == fid] for fid in sorted(families)}
    ordered = list(family_rows)
    for offset in range(0, len(ordered), 20):
        subset = ordered[offset:offset + 20]
        filename = f'family_batches/F{subset[0]}-F{subset[-1]}.md'
        content = f'# Families {subset[0]}-{subset[-1]}\n\nSource: {BASELINE}. Metadata, not proof verification.\n'
        for fid in subset:
            group = family_rows[fid]
            r = group[0]
            content += f"\n## F{fid}. {r['family_title']}\n\nDiscipline: {r['discipline']}\n\n{r['family_description']}\n\n"
            content += '\n'.join(f"- {x['paper_id']}: [{x['title']}]({x['source_url']})" for x in group) + '\n'
        output['our_review/' + filename] = content
        family_lines.append(f'- [{subset[0]}-{subset[-1]}]({filename})\n')
    output['our_review/372_family_batches.md'] = ''.join(family_lines)
    return output

class GitHubAPI:
    def __init__(self, repository: str, token: str):
        if repository != ALLOWED_REPOSITORY:
            raise ValueError('Remote mode is restricted to the authorized research fork.')
        self.base = 'https://api.github.com/repos/' + repository
        self.token = token

    def request(self, path: str, method: str = 'GET', body: dict | None = None):
        data = None if body is None else json.dumps(body, ensure_ascii=False).encode('utf-8')
        headers = {'Accept': 'application/vnd.github+json', 'User-Agent': 'independent-math-reading-ledger', 'X-GitHub-Api-Version': '2022-11-28'}
        if self.token:
            headers['Authorization'] = 'Bearer ' + self.token
        if data is not None:
            headers['Content-Type'] = 'application/json'
        req = urllib.request.Request(self.base + path, data=data, headers=headers, method=method)
        with urllib.request.urlopen(req, timeout=45) as response:
            return json.load(response)

    def file(self, path: str, ref: str) -> bytes:
        result = self.request('/contents/' + urllib.parse.quote(path, safe='/') + '?ref=' + urllib.parse.quote(ref, safe=''))
        if result.get('encoding') != 'base64' or not result.get('content'):
            raise ValueError(f'No complete base64 file content for {path}')
        data = base64.b64decode(result['content'], validate=False)
        if len(data) != result['size'] or git_blob(data) != result['sha']:
            raise ValueError(f'File length or blob hash mismatch: {path}')
        return data

    def publish(self, output: dict[str, str]) -> str:
        head = self.request('/git/ref/heads/main')['object']['sha']
        base_tree = self.request('/git/commits/' + head)['tree']['sha']
        entries = []
        for path, content in output.items():
            if not path.startswith('our_review/') or '..' in Path(path).parts:
                raise ValueError('Refusing write outside generated review directory.')
            entries.append({'path': path, 'mode': '100644', 'type': 'blob', 'content': content})
        tree = self.request('/git/trees', 'POST', {'base_tree': base_tree, 'tree': entries})['sha']
        commit = self.request('/git/commits', 'POST', {'message': 'Generate complete 722-manuscript census; preserve explicit reading status', 'tree': tree, 'parents': [head]})['sha']
        self.request('/git/refs/heads/main', 'PATCH', {'sha': commit, 'force': False})
        return commit

def self_test() -> None:
    fixture = '**001. First family.** Summary. ([Lean](lean/docs/001.md))\n\n&emsp;[First paper](preprints/first/paper.pdf)\n\nFirst abstract.\n</td>\n**002. Second family.** Summary two.\n\n&emsp;[Second paper](preprints/second/paper.pdf)\n\nSecond abstract.\n\n&emsp;[Third paper](preprints/third/main.pdf)\n\nThird abstract.\n'
    rows = parse_catalogue(fixture)
    assert len(rows) == 3
    assert [r['family_id'] for r in rows] == ['001', '002', '002']
    assert rows[0]['formalization_link_listed_for_family'] is True
    assert rows[1]['formalization_link_listed_for_family'] is False
    assert rows[2]['abstract'] == 'Third abstract.'
    assert parse_disciplines('\\cataloguesection{Physics}{1}\n\\resultentry{001}{title}') == {'001': 'Physics'}
    assert git_blob(b'') == 'e69de29bb2d1d6434b8b29ae775ad8c2e48c5391'
    print('Parser self-tests: passed (synthetic fixture, not manuscript verification).')

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default='.')
    parser.add_argument('--remote', action='store_true')
    parser.add_argument('--publish', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.publish and not args.remote:
        raise ValueError('--publish requires --remote.')
    root = Path(args.root)
    if args.remote:
        api = GitHubAPI(os.environ['GITHUB_REPOSITORY'], os.environ.get('GITHUB_TOKEN', ''))
        catalogue = api.file('CONTENTS.md', BASELINE)
        overview = api.file('overview.tex', BASELINE)
        state = json.loads(api.file('our_review/review_state.json', os.environ.get('GITHUB_SHA', 'main')))
    else:
        catalogue = (root / 'CONTENTS.md').read_bytes()
        overview = (root / 'overview.tex').read_bytes()
        state = json.loads((root / 'our_review/review_state.json').read_text(encoding='utf-8'))
    output = render(catalogue, overview, state)
    for path, content in output.items():
        destination = root / path
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(content, encoding='utf-8', newline='\n')
    print(output['our_review/index_audit.json'])
    if args.publish:
        print('Published commit:', api.publish(output))

if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, urllib.error.URLError) as exc:
        print(f'Index build failed: {exc}', file=sys.stderr)
        sys.exit(1)
