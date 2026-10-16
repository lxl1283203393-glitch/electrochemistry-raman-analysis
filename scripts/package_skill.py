"""Create a distributable source skill ZIP, excluding environments and user data."""
from __future__ import annotations
import argparse
import hashlib
import os
from pathlib import Path
import zipfile


def package(root: Path,destination: Path) -> dict:
    root=root.resolve(); destination=destination.resolve()
    if destination.is_relative_to(root): raise ValueError('archive output must be outside the skill source')
    included=[]
    public_root_files={'SKILL.md','README.md','VALIDATION.md','RELEASE_NOTES.md','LICENSE',
        'pyproject.toml','requirements.txt','requirements-tested.txt','.gitignore','.gitattributes'}
    public_directories={'agents','assets','optosky_raman','references','scripts','tests','.github'}
    public_extensions={'.py','.md','.yaml','.yml','.toml','.txt','.svg','.png','.json','.gitignore','.gitattributes'}
    for current,dirs,files in os.walk(root):
        dirs[:]=sorted(d for d in dirs if d not in {'.venv','.git','__pycache__','.pytest_cache','outputs','build','dist'}
            and not d.endswith('.egg-info') and not getattr(Path(current)/d,'is_junction',lambda:False)())
        for name in sorted(files):
            path=Path(current)/name
            relative=path.relative_to(root)
            if len(relative.parts)==1 and name not in public_root_files: continue
            if len(relative.parts)>1 and relative.parts[0] not in public_directories: continue
            if len(relative.parts)>1 and path.suffix.lower() not in public_extensions: continue
            if name in {'results.json','errors.json','verification.json','figures_manifest.json','config_used.yaml'}: continue
            if path.suffix in {'.pyc','.zip'} or name=='MANIFEST.sha256': continue
            included.append(path)
    manifest=''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(root).as_posix()}\n' for p in included)
    (root/'MANIFEST.sha256').write_text(manifest,encoding='utf-8',newline='')
    included.append(root/'MANIFEST.sha256')
    destination.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(destination,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for path in included:
            relative=f'{root.name}/{path.relative_to(root).as_posix()}'
            info=zipfile.ZipInfo(relative,date_time=(1980,1,1,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            archive.writestr(info,path.read_bytes())
    with zipfile.ZipFile(destination) as archive:
        assert archive.testzip() is None
        for line in manifest.splitlines():
            expected,relative=line.split('  ',1)
            assert hashlib.sha256(archive.read(f'{root.name}/{relative}')).hexdigest()==expected
    return {'archive':str(destination),'files':len(included),'bytes':destination.stat().st_size,
        'sha256':hashlib.sha256(destination.read_bytes()).hexdigest()}


if __name__=='__main__':
    import json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    print(json.dumps(package(Path(__file__).resolve().parents[1],args.output),ensure_ascii=False))
