# Installing and Running Leonard

**Leonard** performs high-dimensional Laplacian clustering and Exhaustive Projection Pursuit (EPP) directly from FlowJo workspaces (`.wsp`).

To make it as simple as possible for bench biologists, Leonard is distributed as a **desktop droplet** for both Windows and macOS. You do not need to open a terminal, configure environment variables, or know how to code.

---

## Windows Installation & Usage

There are two options for Windows: the **Setup Installer (Recommended)** or the **Portable Zip**.

### Option A: Setup Installer (Recommended)

1. Download **`Leonard-Setup-x64.exe`** from the latest [GitHub Release](https://github.com/black98fxstc/covariant/releases).
2. Double-click **`Leonard-Setup-x64.exe`** to install.
   * *Note: No administrator privileges are required.* It installs directly to your user folder.
   * If Windows SmartScreen displays *"Windows protected your PC"*, click **More info** &rarr; **Run anyway**.
3. The installer creates a shortcut on your Desktop named **Leonard (Drop .wsp here)** and adds a right-click context menu in File Explorer.

#### How to Analyze an Experiment:
* **Drag-and-Drop (Fastest):** Drag any FlowJo workspace file (`.wsp`) and drop it directly onto the **Leonard** desktop icon.
* **Right-Click Context Menu:** Right-click any `.wsp` file in Windows Explorer and select **Analyze with Leonard**.
* **Double-Click:** Double-click the **Leonard** desktop icon. A file browser window will pop up allowing you to browse to and select your `.wsp` file.

---

### Option B: Portable Folder (No Installation Needed)

1. Download **`Leonard-Windows-Portable.zip`** from the releases page.
2. Extract the `.zip` archive to a folder on your computer (e.g. `Desktop` or `Documents`).
3. To run:
   * Drag your `.wsp` file directly onto **`Leonard.bat`** (or `leonard.exe`).
   * Or double-click **`Leonard.bat`** to open the file selection dialog.

---

## macOS Installation & Usage

On macOS, Leonard is packaged as a standard Mac Application droplet (`Leonard.app`).

### 1. Download & Install

1. Download **`Leonard-macOS.zip`** from the latest [GitHub Release](https://github.com/black98fxstc/covariant/releases).
2. Double-click the `.zip` file to extract it. You will see **`Leonard.app`**.
3. Drag **`Leonard.app`** into your **`Applications`** folder (or onto your Desktop or Dock).

### 2. Security Override (One-Time Setup)

Because this version is distributed outside the Apple App Store, macOS Gatekeeper will protect the app on its first run. You only need to perform this override **once**:

1. Open the **Terminal** app (press `Cmd + Space`, type `Terminal`, and hit `Return`).
2. Paste the following command and hit `Return`:
   ```bash
   xattr -cr /Applications/Leonard.app