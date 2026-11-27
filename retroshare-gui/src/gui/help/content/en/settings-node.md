# Node

The Node page of the [Preferences](settings.md) shows the identifiers of your
node and of your profile, and lets you export your profile and your friends
list. Nothing on this page is a setting: it is information and three buttons.

Your RetroShare ID, the text you give to a friend, is on the
[Home](home.md) page.

## Layout

The page is split into three blocks.

**Public Information**

- **Name:** the name of your profile.
- **Location:** the name of this node.
- **Location ID:** the identifier of this node. Each node of a profile has
  its own.
- **Software Version:** shown as **RetroShare version:** and **libretroshare
  version:**, the version of the program and of its core library.
- **Online since:** the date and time RetroShare was started.

**Profile**

- **PGP Id:** the identifier of your profile key.
- **PGP fingerprint:** the fingerprint of your profile key. Compare it with a
  friend over another channel to check that a profile really is yours.
- **Friend nodes:** the number of profiles you have accepted as friends.
- **Export Profile:** the **Export profile** button.
- **Export Friendslist:** the **Export friends list** button.

**Other Information**

- **Profile path:** the directory where the data of this node is stored.
- **Statistics:** the **Show statistics** button.

## Buttons

- **Export profile**: saves your profile key pair to an `.asc` file. The file
  is encrypted: you need the password of your profile to use it. Import it on
  another computer to create a new node for the same profile; the friends you
  add on the new node recognise it as yours. RetroShare confirms with
  **Your identity was successfully saved**, or reports that an error
  occurred.
- **Export friends list**: saves the list of your friends and of their groups
  to an `.xml` file that you can import on another node. RetroShare reminds
  you that this file is not encrypted.
- **Show statistics**: opens the [Statistics](statistics.md) window.

## How to

### Use your profile on a second computer

1. Click **Export profile** and choose where to save the `.asc` file.
2. Copy the file to the other computer.
3. On the other computer, start RetroShare and import the file when you
   create the profile. See [Profiles and nodes](start-profile.md).
4. Optionally, click **Export friends list** here and import the `.xml` file
   on the new node to get the same friends.

## See also

- [Home](home.md)
- [Profiles and nodes](start-profile.md)
- [Network](settings-network.md)
