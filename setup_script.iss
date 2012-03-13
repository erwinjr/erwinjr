; Inno Setup Script for version 5.4.2

[Setup]
AppName=ErwinJr
AppVersion=2.1
DefaultDirName={pf}\ErwinJr
DefaultGroupName=ErwinJr
UninstallDisplayIcon={app}\EJico.ico
Compression=lzma2
SolidCompression=yes
OutputDir=c:\erwinjr
LicenseFile=license.txt
ChangesAssociations=yes

[Files]
Source: "dist/*"; DestDir: "{app}"
Source: "dist/compiler/*"; DestDir: "{app}/compiler"
Source: "dist/ctypes/*"; DestDir: "{app}/ctypes"
Source: "dist/dateutil/*"; DestDir: "{app}/dateutil"
Source: "dist/dateutil/zoneinfo/*"; DestDir: "{app}/dateutil/zoneinfo"
Source: "dist/distutils/*"; DestDir: "{app}/distutils"
Source: "dist/email/*"; DestDir: "{app}/email"
Source: "dist/email/mime/*"; DestDir: "{app}/email/mime"
Source: "dist/encodings/*"; DestDir: "{app}/encodings"
Source: "dist/examples/*"; DestDir: "{app}/examples"
Source: "dist/images/*"; DestDir: "{app}/images"
Source: "dist/logging/*"; DestDir: "{app}/logging"
Source: "dist/matplotlib/*"; DestDir: "{app}/matplotlib"
Source: "dist/matplotlib/backends/*"; DestDir: "{app}/matplotlib/backends"
Source: "dist/matplotlib/backends/qt4_editor/*"; DestDir: "{app}/matplotlib/backends/qt4_editor"
Source: "dist/matplotlib/delaunay/*"; DestDir: "{app}/matplotlib/delaunay"
Source: "dist/matplotlib/projections/*"; DestDir: "{app}/matplotlib/projections"
Source: "dist/matplotlib/testing/*"; DestDir: "{app}/matplotlib/testing"
Source: "dist/matplotlib/tri/*"; DestDir: "{app}/matplotlib/tri"
Source: "dist/mpl_toolkits/*"; DestDir: "{app}/mpl_toolkits"
Source: "dist/mpl_toolkits/mplot3d/*"; DestDir: "{app}/mpl_toolkits/mplot3d"
Source: "dist/mpl-data/*"; DestDir: "{app}/mpl-data"
Source: "dist/mpl-data/fonts/afm/*"; DestDir: "{app}/mpl-data/fonts/afm"
Source: "dist/mpl-data/fonts/pdfcorefonts/*"; DestDir: "{app}/mpl-data/fonts/pdfcorefonts"
Source: "dist/mpl-data/fonts/ttf/*"; DestDir: "{app}/mpl-data/fonts/ttf"
Source: "dist/mpl-data/images/*"; DestDir: "{app}/mpl-data/images"
Source: "dist/multiprocessing/*"; DestDir: "{app}/multiprocessing"
Source: "dist/multiprocessing/dummy/*"; DestDir: "{app}/multiprocessing/dummy"
Source: "dist/nose/*"; DestDir: "{app}/nose"
Source: "dist/nose/ext/*"; DestDir: "{app}/nose/ext"
Source: "dist/nose/plugins/*"; DestDir: "{app}/nose/plugins"
Source: "dist/numpy/*"; DestDir: "{app}/numpy"
Source: "dist/numpy/compat/*"; DestDir: "{app}/numpy/compat"
Source: "dist/numpy/core/*"; DestDir: "{app}/numpy/core"
Source: "dist/numpy/fft/*"; DestDir: "{app}/numpy/fft"
Source: "dist/numpy/lib/*"; DestDir: "{app}/numpy/lib"
Source: "dist/numpy/linalg/*"; DestDir: "{app}/numpy/linalg"
Source: "dist/numpy/ma/*"; DestDir: "{app}/numpy/ma"
Source: "dist/numpy/matrixlib/*"; DestDir: "{app}/numpy/matrixlib"
Source: "dist/numpy/numarray/*"; DestDir: "{app}/numpy/numarray"
Source: "dist/numpy/oldnumeric/*"; DestDir: "{app}/numpy/oldnumeric"
Source: "dist/numpy/polynomial/*"; DestDir: "{app}/numpy/polynomial"
Source: "dist/numpy/random/*"; DestDir: "{app}/numpy/random"
Source: "dist/numpy/testing/*"; DestDir: "{app}/numpy/testing"
Source: "dist/psyco/*"; DestDir: "{app}/psyco"
Source: "dist/PyQt4/*"; DestDir: "{app}/PyQt4"
Source: "dist/PyQt4/Qwt5/*"; DestDir: "{app}/PyQt4/Qwt5"
Source: "dist/pyreadline/*"; DestDir: "{app}/pyreadline"
Source: "dist/pyreadline/clipboard/*"; DestDir: "{app}/pyreadline/clipboard"
Source: "dist/pyreadline/console/*"; DestDir: "{app}/pyreadline/console"
Source: "dist/pyreadline/keysyms/*"; DestDir: "{app}/pyreadline/keysyms"
Source: "dist/pyreadline/lineeditor/*"; DestDir: "{app}/pyreadline/lineeditor"
Source: "dist/pyreadline/modes/*"; DestDir: "{app}/pyreadline/modes"
Source: "dist/pytz/*"; DestDir: "{app}/pytz"
Source: "dist/scipy/*"; DestDir: "{app}/scipy"
Source: "dist/scipy/interpolate/*"; DestDir: "{app}/scipy/interpolate"
Source: "dist/scipy/linalg/*"; DestDir: "{app}/scipy/linalg"
Source: "dist/scipy/misc/*"; DestDir: "{app}/scipy/misc"
Source: "dist/scipy/sparse/*"; DestDir: "{app}/scipy/sparse"
Source: "dist/scipy/spatial/*"; DestDir: "{app}/scipy/spatial"
Source: "dist/scipy/special/*"; DestDir: "{app}/scipy/special"
Source: "dist/scipy/sparse/sparsetools/*"; DestDir: "{app}/scipy/sparse/sparsetools"
Source: "dist/src/*"; DestDir: "{app}/src"
;Source: "dist/src/images/*"; DestDir: "{app}/src/images"
;Source: "dist/wx/*"; DestDir: "{app}/wx"
Source: "dist/xml/*"; DestDir: "{app}/xml"
Source: "dist/xml/parsers/*"; DestDir: "{app}/xml/parsers"
Source: "dist/xml/sax/*"; DestDir: "{app}/xml/sax"

[Icons]
Name: "{commonprograms}\ErwinJr\ErwinJr"; Filename: "{app}\erwinjr.exe"; IconFilename: "{app}\EJico.ico"; WorkingDir: "{app}"
Name: "{commonprograms}\ErwinJr\Uninstall ErwinJr"; Filename: "{uninstallexe}"
Name: "{commondesktop}\ErwinJr"; Filename: "{app}\erwinjr.exe"; IconFilename: "{app}\EJico.ico"; WorkingDir: "{app}"

[Registry]
Root: HKCR; Subkey: ".qcl"; ValueType: string; ValueName: ""; ValueData: "QCL.Document"; Flags: uninsdeletevalue
Root: HKCR; Subkey: "QCL.Document"; ValueType: string; ValueName: ""; ValueData: "ErwinJr qcl file"; Flags: uninsdeletekey
Root: HKCR; Subkey: "QCL.Document\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}/EJico.ico"
Root: HKCR; Subkey: "QCL.Document\shell\Open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\erwinjr.exe"" ""%1"""
Root: HKCU; Subkey: "Software\JPL\ErwinJr"; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\JPL"; Flags: uninsdeletekeyifempty
Root: HKLM; Subkey: "Software\JPL\ErwinJr"; Flags: uninsdeletekey
Root: HKLM; Subkey: "Software\JPL"; Flags: uninsdeletekeyifempty
;Root: HKLM; Subkey: "Software\Wow6432Node\JPL\ErwinJr"; Flags: uninsdeletekey
;Root: HKLM; Subkey: "Software\Wow6432Node\JPL"; Flags: uninsdeletekeyifempty

[Run]
Filename: "{app}\erwinjr.exe"; Description: "Launch application"; Flags: postinstall nowait skipifsilent
