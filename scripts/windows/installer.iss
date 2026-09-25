#define MyAppName "Leonard"
#define MyAppVersion "0.1.0"
#define MyAppPublisher "Covariant"
#define MyAppExeName "leonard.exe"
#define MyAppCmdName "Leonard.cmd"

[Setup]
AppId={{D37F2936-7C1C-4B14-B1F7-418E66275811}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={localappdata}\Programs\{#MyAppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir=..\..\dist
OutputBaseFilename=Leonard-Setup-x64
Compression=lzma2
SolidCompression=yes
WizardStyle=modern

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Files]
Source: "..\..\build\Release\{#MyAppExeName}"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\..\build\Release\*.dll"; DestDir: "{app}"; Flags: ignoreversion skipifsourcedoesntexist
Source: "{#MyAppCmdName}"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{autoprograms}\{#MyAppName}"; Filename: "{app}\{#MyAppCmdName}"; IconFilename: "{app}\{#MyAppExeName}"
Name: "{userdesktop}\{#MyAppName} (Drop .wsp here)"; Filename: "{app}\{#MyAppCmdName}"; IconFilename: "{app}\{#MyAppExeName}"

[Registry]
; File Explorer Right-Click Context Menu for .wsp files
Root: HKCU; Subkey: "Software\Classes\SystemFileAssociations\.wsp\shell\Leonard"; ValueType: string; ValueData: "Analyze with Leonard"; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\Classes\SystemFileAssociations\.wsp\shell\Leonard\command"; ValueType: string; ValueData: """{app}\{#MyAppCmdName}"" ""%1"""

; Fallback for FlowJo ProgID association if registered
Root: HKCU; Subkey: "Software\Classes\FlowJo.Workspace\shell\Leonard"; ValueType: string; ValueData: "Analyze with Leonard"; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\Classes\FlowJo.Workspace\shell\Leonard\command"; ValueType: string; ValueData: """{app}\{#MyAppCmdName}"" ""%1"""
