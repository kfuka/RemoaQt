# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['richgui.py'],
             pathex=['D:\\GoogleDrive\\Kanagawa\\20200709\\Final\\Remoa_Qt'],
             binaries=[],
             datas=[],
             hiddenimports=['statsmodels.tsa.statespace._kalman_filter', 'statsmodels.tsa.statespace._kalman_smoother', 'statsmodels.tsa.statespace._representation', 'statsmodels.tsa.statespace._simulation_smoother', 'statsmodels.tsa.statespace._statespace', 'statsmodels.tsa.statespace._tools', 'statsmodels.tsa.statespace._filters._conventional', 'statsmodels.tsa.statespace._filters._inversions', 'statsmodels.tsa.statespace._filters._univariate', 'statsmodels.tsa.statespace._smoothers._alternative', 'statsmodels.tsa.statespace._smoothers._classical', 'statsmodels.tsa.statespace._smoothers._conventional', 'statsmodels.tsa.statespace._smoothers._univariate', 'statsmodels.tsa.statespace._filters.univariate_diffuse', 'statsmodels.tsa.statespace._smoothers._univariate_diffuse'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='richgui',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True )
