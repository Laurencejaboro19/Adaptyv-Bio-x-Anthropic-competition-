# Block 4: download a finished run from the Modal volume and unpack it. Set $tag to the tag you launched.
$dir = "$env:USERPROFILE\Downloads\tfr1_rfd4"
$tag = "my_first_run"
Set-Location $dir
$env:PYTHONUTF8="1"; $env:PYTHONIOENCODING="utf-8"
python -m modal run rfd4_app.py::fetch --tag $tag --outdir "$dir\out"
tar -xzf "$dir\out\$tag.tar.gz" -C "$dir\out"
Get-ChildItem "$dir\out\$tag" -Recurse -File | Select-Object -ExpandProperty FullName
