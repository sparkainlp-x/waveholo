# Third-party data

The code in this repository is AGPL-3.0-only (see `LICENSE`). The
following third-party data is included or derived from, under its own
licence.

## pyDHM fly-head hologram (MIT)

- File: `real_data/fly-20x.bmp`, copied unmodified from
  https://github.com/catrujilla/pyDHM (`data/compensation samples/`).
- Reconstructions in `real_data/fly_*.png` and `results/s10_*.png` are
  derived from it.
- Licence text, as published by pyDHM:

```
MIT License

Copyright (c) 2022 Ana Doblas and Carlos Alejandro Trujillo (EAFIT)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Cuenat et al. DHM dataset (CC-BY-4.0)

- Source: Cuenat et al., "Fast Autofocusing using Tiny Networks for
  Digital Holographic Microscopy" (arXiv:2203.07772), dataset Zenodo
  record 6337535,
  https://doi.org/10.5281/zenodo.6337535, licensed CC-BY-4.0
  (https://creativecommons.org/licenses/by/4.0/).
- No raw holograms from the dataset are redistributed here. The
  figures and tables in `real_data/dhm_tiny/` and `real_data/dhm_v013/`
  are derived from it: holograms were cropped, reconstructed,
  numerically refocused and summarised by this package.

## Not redistributed

- `ODP-DLHM-Database.csv` (the HQ-DLHM-OD-Database metadata index,
  https://doi.org/10.17605/OSF.IO/SYHT9) is referred to in
  `real_data/REAL_DATA.md` but is not included in this public
  repository, because its OSF licence is listed as "Other" with no
  terms stated. Download it from the OSF page. No script uses it.
