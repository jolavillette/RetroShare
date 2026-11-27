# General

The General page of the [Preferences](settings.md) sets how RetroShare starts,
what happens when you close its window and when your status becomes inactive.

The language and the look of the program are in
[Appearance](settings-appearance.md).

## Layout

The page holds three groups: **Startup**, **Misc** and **Idle**.

## Buttons

### Startup

- **Start RetroShare when my system starts**: only on Windows. RetroShare is
  launched when you open your Windows session.
- **Start minimized**: RetroShare starts without showing its main window. Use
  the tray icon to show it. Off by default.
- **Start minimized on system start**: only on Windows. When RetroShare is
  launched by the system start, it starts minimized.

### Misc

- **Do not show the Quit RetroShare MessageBox**: RetroShare quits at once
  when you choose **Quit**. When the box is not checked, RetroShare first asks
  **Do you really want to exit RetroShare ?**. Off by default.
- **Minimize to Tray Icon**: closing the main window hides it and RetroShare
  keeps running; the tray icon brings the window back. When the box is not
  checked, closing the window quits RetroShare. If your system has no tray
  that RetroShare can use, the option is switched off and its tooltip says
  **No Qt-compatible system tray was found on this system.**
- **Auto Login**: RetroShare starts your node without asking for the password.
  The tooltip warns that this is discouraged for security reasons. The box is
  greyed, with the tooltip **Your RetroShare build has auto-login disabled.**,
  when RetroShare was built without this feature.
- **Register retroshare:// as URL protocol**: only on Windows. Your system
  then opens `retroshare://` links, for instance in a web page or an e-mail,
  with RetroShare. See [RetroShare links](links.md). The box is greyed when
  your account is not allowed to change this; the icon next to it then has the
  tooltip **You don't have sufficient rights. Run RetroShare as Admin to
  change this setting.** If the change fails, a message says **Could not add
  retroshare:// as protocol.** or **Could not remove retroshare:// protocol.**
- **Accept operating systems calls**: RetroShare accepts the requests of your
  operating system to open RetroShare collection files and download links.
  This is what makes a click on a `retroshare://` link outside RetroShare
  reach the running program. On by default.
- **!!!The RetroShare's desktop file is missing or wrong!!!**: not on Windows.
  This warning is shown when the `retroshare://` protocol is not registered on
  your system. Its tooltip says to install RetroShare with a package
  installer to get `/usr/share/applications/retroshare.desktop`.

### Idle

- **Idle Time**: the time without keyboard or mouse activity after which your
  status is set to "Inactive". In seconds, from 0 to 1000, 300 by default. The
  value 0 is shown as **Never** and disables the feature.

## See also

- [Main window, tray icon and status bar](mainwindow.md)
- [Command line options](commandline.md)
