$work = "$env:USERPROFILE\Downloads\glp1r_refold"
$cifs = "$work\cifs"
$arch = "$work\archives"
New-Item -ItemType Directory -Force -Path $cifs, $arch | Out-Null
$ids = @'
bridge_01 sab_pred_lQ48TIB9pnLtI7deIZv3
bridge_02 sab_pred_zx1aLspbbliYdIchlT3N
bridge_03 sab_pred_W2RhFePW5MjHF1xzKm6a
bridge_04 sab_pred_U97HrKCXbgGaUXD9ANyW
bridge_05 sab_pred_qUw2B94kwqfHSu1K23bK
bridge_06 sab_pred_LAyF1Q2Z9n130Dv4gSHE
bridge_07 sab_pred_9GdYFJM6th8t8XdKoC2l
bridge_08 sab_pred_nVzIGM8BEpoNDn3xcvWe
bridge_09 sab_pred_TpMsScJHQeUK9ZWLjUfm
bridge_10 sab_pred_3RM7NjKvcV3yGJHoGmWW
bridge_11 sab_pred_WhMMHaPt3MEZEKzcyGat
bridge_12 sab_pred_2qlsxA9kF7p4hyLUlUWw
bridge_13 sab_pred_vDyo4ytaVbeK8dNwaVMl
bridge_14 sab_pred_ftWcorJR0hXrl1SOR1D5
bridge_15 sab_pred_3aaT2lZyimuyHB4CcP5J
bridge_16 sab_pred_1Aj08MQLcfY6rkXYm2Pv
bridge_17 sab_pred_UCesPjoUfAff3RNmE8Wx
bridge_18 sab_pred_hgKsj79VZJoFuV326hEp
bridge_19 sab_pred_LSwHIN4hqMh9YemnfCWn
bridge_20 sab_pred_ezB55TNVHFhrKxfnCIlN
bridge_21 sab_pred_GQxQJf7ZA6fSo43V8LUA
bridge_22 sab_pred_fttflr0eucTWpzT36aAv
bridge_23 sab_pred_489DGITYM6woV31hZoDT
bridge_24 sab_pred_E00jRqbaAaOtObjaeCKm
bridge_25 sab_pred_MpfcyCK6R7nfkxHIgu6x
bridge_26 sab_pred_T7RtOH2js7vJPHu0c5th
bridge_27 sab_pred_Y8rrGtc5HgTvRy6ft10Y
bridge_28 sab_pred_B37ztZJK8gv07RDgPOJI
bridge_29 sab_pred_Jigi7aXsyrm6GS6m5N71
bridge_30 sab_pred_n6SpBRnOW4Ih7mFzWMNK
bridge_31 sab_pred_5TVd6acPc1OOXRRMvDqY
bridge_32 sab_pred_U0fV58E9IwUEmnE8ZNt1
bridge_33 sab_pred_gsUq1TT2jDTDze9S5Goj
bridge_34 sab_pred_vWotpIekILrylQYGeLkI
bridge_35 sab_pred_SZEoO1VZCGkwh3Q7MDtd
bridge_36 sab_pred_AF720bn52mYq20h0clnw
bridge_37 sab_pred_jCGqC9OVSP4qcDEoNjR3
bridge_38 sab_pred_cuOfXj6x5sne3gDHFlDJ
bridge_39 sab_pred_u7SzLUhVc9GXx7xEPBPn
bridge_40 sab_pred_20mxISTOidYLOMEDssfL
bridge_41 sab_pred_FBBJjURVlZmXD0RjKDEF
bridge_42 sab_pred_M7XtWEpFCLkkCCym6Qdo
bridge_43 sab_pred_oP3NL9TjtBdFPOGbo5tO
bridge_44 sab_pred_3dzFvfc8SsWtpjfmvk8r
bridge_45 sab_pred_i0meMThozb87ZaG4nYx8
bridge_46 sab_pred_koUECGCO0KCs0rEY0dGu
bridge_47 sab_pred_VwOOV00fyzsH5BMuCJWu
bridge_48 sab_pred_i8a6BLWxmFcwgRcKCGWN
bridge_49 sab_pred_0zpZl4k2BbIFjYSLQ0Dh
bridge_50 sab_pred_LclfbKQsgi5tLivqMIa1
deep_01 sab_pred_yUdrOyQxh6mlBxKvlUnN
deep_02 sab_pred_6mxNKk3OhhzKE5dwI2bo
deep_03 sab_pred_pCcFayj6iTPcjcjlupxx
deep_04 sab_pred_vuTGasl927TYJUTudwTB
deep_05 sab_pred_Xt1X77LuWz3v9wne0DHX
deep_06 sab_pred_09RMAJliBkbDgIKVOSpz
deep_07 sab_pred_gYi9bA04ottXBsjfiJju
deep_08 sab_pred_RCfwZFCZbwOXd8eIdbEr
deep_09 sab_pred_eYtAzrSPVrFpgTem2VnC
deep_10 sab_pred_hd1vhyQnDVGxJCE9PXVz
deep_11 sab_pred_xX7cUqU4HxDVFwbnMFUi
deep_12 sab_pred_qJCosSGGskWkDki4lx91
deep_13 sab_pred_WEAiwzwzdyGHryPJfYcN
deep_14 sab_pred_MieShfxFWLZrYD0nj7Kg
deep_15 sab_pred_1y1jsnX0bRNRrRgbyd6S
deep_16 sab_pred_e1KcjxgqhoJK8C36qp3u
deep_17 sab_pred_nKfRFNgzm29puXDCTes5
deep_18 sab_pred_u8p4QZkOicY8qdhsG7UD
deep_19 sab_pred_TuHPY8sCMPn1d9P4oZJ2
deep_20 sab_pred_YMDAmnVoFznfkK5HSPQV
deep_21 sab_pred_sOeHw9XV2GBtDjrwEcLY
deep_22 sab_pred_3W2ohQqCrrHKXFa7YdKi
deep_23 sab_pred_4wWZMMP0rP5HgurDZvcL
deep_24 sab_pred_I0cIko2NKaX9mxt1L7sk
deep_25 sab_pred_QZENreMfyUkln5esRM3u
deep_26 sab_pred_6lAPaDyBCcLIFx6Pa9LC
deep_27 sab_pred_iCOu6gQMeOcjwRG7LHFb
deep_28 sab_pred_oeXENgtRRXNzu1HEpwGg
deep_29 sab_pred_gNWpMng7SAxxG6r7rMCM
deep_30 sab_pred_wzPWlK3K17gEGbqwLXJg
deep_31 sab_pred_B5uNuVKGz1XNDEy1P2dK
deep_32 sab_pred_cEf3qvVYWy5Bl0ZZgISS
deep_33 sab_pred_iU1qGHbkZC6vqcGM4Cji
deep_34 sab_pred_pZkgB2jLnK3pqsRloTmk
deep_35 sab_pred_nzWi5p4x3366hqNgRGRH
deep_36 sab_pred_oNKzK6Gu2FOLWOaSuqzW
deep_37 sab_pred_qhdW1oisJgC09353aq8D
deep_38 sab_pred_n3fGPqUtuAI3WoaxYQ5v
deep_39 sab_pred_QThYnrPFO7IhnHymIR2O
deep_40 sab_pred_COIh8JUFWe4Gx9GuXivA
deep_41 sab_pred_eFhnEInbHt9nRCU0y39P
deep_42 sab_pred_moyWvx4PgquFAMepVqNn
deep_43 sab_pred_r92VxFFQzEOv85bfBHxH
deep_44 sab_pred_HxBdGeeAOvupUufUr3HJ
deep_45 sab_pred_deK7RGFKuwno3a6qmpwK
deep_46 sab_pred_feWS8Q3RhRjYzZdocRPT
deep_47 sab_pred_XM87aRsueGYPWEwbsJVU
deep_48 sab_pred_1tcrGUAZdMg5bnOTLagJ
deep_49 sab_pred_21yKvJG69vZ3m7LmjXWp
deep_50 sab_pred_xrOXhx6zPCjjEL8mAR4H
'@
$env:BOLTZ_API_KEY = $null
$sec = Read-Host "Paste your Boltz API key (right-click to paste, then Enter)" -AsSecureString
$plain = [System.Net.NetworkCredential]::new("", $sec).Password
$env:BOLTZ_API_KEY = ($plain -replace '[^\x21-\x7E]', '')
Write-Host ("Key characters received: " + $env:BOLTZ_API_KEY.Length)
if ($env:BOLTZ_API_KEY.Length -lt 10) { Write-Host "Key looks too short - paste did not work. Stop here and try again."; return }
$hdr = @{ "x-api-key" = $env:BOLTZ_API_KEY }
$test = Invoke-RestMethod -Headers $hdr -Uri "https://api.boltz.bio/compute/v1/predictions/structure-and-binding/sab_pred_lQ48TIB9pnLtI7deIZv3"
Write-Host ("Test call status: " + $test.status)
$fail = @()
foreach ($line in ($ids -split "`r?`n")) {
  if (-not $line.Trim()) { continue }
  $p = $line.Trim() -split "\s+"
  $name = $p[0]; $id = $p[1]
  try {
    $r = Invoke-RestMethod -Headers $hdr -Uri "https://api.boltz.bio/compute/v1/predictions/structure-and-binding/$id"
    Invoke-WebRequest -Uri $r.output.best_sample.structure.url -OutFile "$cifs\$name.cif"
    Invoke-WebRequest -Uri $r.output.archive.url -OutFile "$arch\$name.tar.gz"
    Write-Host "ok $name"
  } catch {
    Write-Host "FAILED $name : $($_.Exception.Message)"
    $fail += $name
  }
}
Write-Host ("cifs downloaded: " + (Get-ChildItem $cifs -Filter *.cif).Count + " / 100; failed: " + $fail.Count)
$ex = "$work\extracted\bridge_01"
New-Item -ItemType Directory -Force -Path $ex | Out-Null
tar -xzf "$arch\bridge_01.tar.gz" -C $ex
Write-Host "Files inside one archive:"
Get-ChildItem $ex -Recurse -File | ForEach-Object { $_.FullName.Substring($ex.Length) }
@'
import sys, glob, os, csv
import gemmi, numpy as np
OFFSET = 23
BRIDGE_HOT = {32,39,69,88,91,121,197,201,205,210,211,214,215,221,298,299,300,306}
DEEP_HOT = {298,299,300,306,230,233,234,309,310}
COLD = {63,82,115,146,207,208,279,282,283,286,290,303,304,307,311,371,373,374,376,378,381,382}
MIN_HOT_FRAC = 0.30
MAX_COLD = 2
MIN_BINDER_CONTACTS = 10
folder = sys.argv[1]
rows = []
for p in sorted(glob.glob(os.path.join(folder, '*.cif'))):
    st = gemmi.read_structure(p); st.setup_entities()
    m = st[0]
    rec, bnd = m['A'], m['B']
    names = {r.seqid.num + OFFSET: r.name for r in rec}
    if names.get(29) != 'THR' or names.get(298) is None:
        sys.exit('Numbering check failed for ' + p + ' (residue 29 should be THR): ' + str(names.get(29)))
    rec_xyz = [(r.seqid.num + OFFSET, np.array([a.pos.tolist() for a in r])) for r in rec]
    ba = [(i + 1, np.array([a.pos.tolist() for a in r])) for i, r in enumerate(bnd)]
    contact_rec, contact_bnd = set(), set()
    for bi, bx in ba:
        for u, rx in rec_xyz:
            dmin = np.linalg.norm(bx[:, None, :] - rx[None, :, :], axis=2).min()
            if dmin <= 4.5:
                contact_rec.add(u); contact_bnd.add(bi)
    L = len(ba)
    hot = BRIDGE_HOT if os.path.basename(p).startswith('bridge') else DEEP_HOT
    hit = len(contact_rec & hot); cold = len(contact_rec & COLD)
    cterm = any(i in contact_bnd for i in range(L - 2, L + 1))
    ok = (hit / len(hot) >= MIN_HOT_FRAC and cold <= MAX_COLD
          and len(contact_bnd) >= MIN_BINDER_CONTACTS and not cterm)
    rows.append(dict(design=os.path.basename(p)[:-4], binder_len=L,
        n_rec_contacts=len(contact_rec), n_binder_contacts=len(contact_bnd),
        hotspots_hit=hit, hotspots_total=len(hot), hotspot_frac=round(hit / len(hot), 2),
        coldspots_hit=cold, cterm_in_interface=cterm, PASS=ok,
        contacted_residues=' '.join(str(x) for x in sorted(contact_rec))))
if not rows:
    sys.exit('No .cif files found in ' + folder)
out = os.path.join(folder, '..', 'contact_check.csv')
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print('wrote', os.path.abspath(out), len(rows), 'designs;', sum(r['PASS'] for r in rows), 'pass')
for c in ('bridge', 'deep'):
    rr = [r for r in rows if r['design'].startswith(c)]
    print(c, len(rr), 'pass', sum(r['PASS'] for r in rr))
'@ | Set-Content -Encoding utf8 "$work\contacts.py"
python -m pip install gemmi numpy
python "$work\contacts.py" $cifs
$env:BOLTZ_API_KEY = $null
