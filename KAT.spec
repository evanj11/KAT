# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_submodules

hiddenimports = (
    collect_submodules('scipy') +
    collect_submodules('sympy') +
    collect_submodules('pandas') +
    collect_submodules('matplotlib') +
    collect_submodules('scikit-learn')
)

a = Analysis(
    ['main.py'],
    pathex=[],
    binaries=[],
    datas=[('scripts', './scripts'), ('assets', './assets'), ('run_script.py', 'scripts'), ('.imgs', './.imgs')],
    hiddenimports=hiddenimports,
    #hiddenimports=['scipy', 'scipy.optimize',
    #'scipy.special',
    #'sympy',
    #'matplotlib',
    #'mpl_toolkits.axes_grid1',
    #'mpl_toolkits.axes_grid1.inset_locator',
    #'pandas',
    #],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'PyQt5',
        'PyQt6',
    ],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='KAT',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['assets/kat.icns'],
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='KAT',
)
app = BUNDLE(
    coll,
    name='KAT.app',
    icon='assets/kat.icns',
    bundle_identifier=None,
)
