# Boards

A board is a place where people share links, pictures and short texts, and
vote on them. Each post has a score made of the votes of its readers, and the
posts can be sorted by score, so the content people like rises to the top.
Every post can be commented.

Boards spread from friend to friend, like [channels](channels.md) and
[forums](forums.md). Anybody who is subscribed to a board can post in it.

## Layout

- The title bar shows **Boards** and the **?** help button.
- On the left, the list of boards with its toolbar.
- On the right, the board that is selected in the list. Several boards can be
  open in tabs.

### The list of boards

The list works like the list of channels; see [Channels](channels.md) for the
columns, the filter field with its **Title** and **Description** options, and
the tooltips. The differences are:

- the categories are **My Boards**, **Subscribed Boards**, **Popular Boards**
  and **Other Boards**;
- the button of the toolbar is **Create Board**;
- there is no field to search the whole network.

### An open board

Two tabs are always present, and one more is added for each post whose
comments you open:

- **Details**: the logo and name of the board with the subscribe button, and
  under **Board Details**: **Popularity:** (the number of your friend nodes
  that are subscribed), **Items (at friends):** (the largest number of posts,
  comments and votes one of your friends has), **Created**,
  **Administrator:**, **Last activity:**, **Distribution:**, **Sync
  period:** (when you are subscribed) and the description.
- **Posts**: a toolbar and the posts, 25 per page.
- A tab named after a post: the comments on that post. Close it with its
  cross.

### A post

Each post shows, from left to right:

- the vote buttons with the score between them;
- the picture of the post, or a default picture;
- the title. When the post carries a link, the title is underlined: click it
  to open the link. Hover it to see the address;
- **Posted by**, the author and the date of the post;
- the buttons described below, and the **New** mark on a post that arrived
  since you last looked.

The score is the number of positive votes minus the number of negative
votes.

A post whose author you banned shows a red notice instead of its title, and
its buttons are disabled.

Two presentations exist: the compact view, with a small picture and one post
per line, and the card view, with a large picture.

## Buttons

### Above the list

- **?** (help button): shows a short introduction to boards.
- **Create Board**: opens the **Create New Board** dialog.

### Details tab

- **Subscribe**: subscribes to the board. Once subscribed, the button reads
  **Subscribed** followed by the popularity, and opens a menu with
  **Unsubscribe**.

### Posts tab

- **Create Post**: opens the **Create a new Post** window. Enabled only when
  you are subscribed.
- Sort box: **New** shows the most recent posts first, **Top** the posts with
  the best score, **Hot** the posts with a good score that are still recent:
  the weight of the score drops as the post gets older.
- **Search**: keeps the posts whose title contains all the words you type.
- View button (**Click to switch to card view** / **Click to switch to
  compact view**): changes the presentation of the posts.
- **Previous**, **Next**: move to the previous or next page of 25 posts. The
  label between them gives the range that is displayed.
- Identity box (tooltip **Default identity used when voting**): the identity
  that signs your votes and is proposed for your comments.

### On each post

- **Vote up**, **Vote down**: adds or removes one point. The vote is signed by
  the identity chosen in the identity box. Once you have voted on a post, its
  two vote buttons are disabled.
- **Comments**: opens the comments of the post in a new tab. The button shows
  the number of comments.
- **Expand**: shows or hides the text of the post under the title. Shown only
  when the post has a text.
- **Toggle Message Read Status**: marks the post as read or unread.
- **Click to view picture** (the picture itself): opens the picture at full
  size in a viewer.

Voting, opening the picture, opening the comments and expanding the text also
mark the post as read.

### Comments tab

- **Comment**: opens the **Make Comment** dialog to write a comment on the
  post.
- Sort box: **Hot**, **New** or **Top** orders the comments.
- **Voter ID:** the identity that signs your votes on comments.

## Context menus

### A board in the list

The entries are the same as for a channel and are described in
[Channels](channels.md):

- **Open in new tab**
- **Subscribe** / **Unsubscribe**
- **Create Board**
- **Show Details**
- **Edit Details**: enabled only for the boards you administrate.
- **Synchronise posts of last...** and **Store posts for at most...**: enabled
  only when you are subscribed.
- **Copy RetroShare Link**
- **Mark all as read**, **Mark all as unread**: enabled only when you are
  subscribed.
- **Tree View Options**: **Show Header**, **Show column …** and, when the
  header is hidden, **Sort by column …**.

A board has no **Share publish permissions...** entry and no download
entries.

### A post

- **Copy RetroShare Link**: copies a link to the post. A message confirms the
  copy.
- **Copy http Link**: copies the web address carried by the post. Only when
  the post has an `http` or `https` link.
- **Show author in People tab**: shows the author in the
  [People](people.md) page.

### A comment

- **Reply to Comment**: opens the **Reply to Comment** dialog.
- **Submit Comment**: writes a new comment on the post.
- **Copy Comment**: copies the text of the comment.
- **Vote Up**, **Vote Down**: votes on the comment.
- **Show Author**: shows the author of the comment in the
  [People](people.md) page.

## Dialogs

### Create New Board

- The picture at the top left is the logo. Click it to choose an image file.
- **Name**: the name of the board. It is required.
- **Owner:** the identity of yours that is shown as the author of the board.
- **Description**: what the board is about.
- **Message Distribution**: who can see the board.
  - **Public**: everybody.
  - **Restricted to Circle**: only the members of the circle you choose. See
    [People](people.md) for circles.
  - **Restricted node group**: only your friend nodes that belong to the
    group you choose.
- **Create**: creates the board. It appears under **My Boards**.

**Edit Details** opens the same dialog, titled **Edit Board**, with the button
**Update Board**. **Show Details** shows the same information read only, with
the ID, the creation date, the number of posts, the last post and the
popularity.

### Create a new Post

- **Title**: the title of the post. It is required.
- **Post as**: the identity that signs the post. You need an identity to post.
- Three buttons choose the page that is displayed. A post can use the three
  pages at once: a text, a picture and a link.
  - **Post**: an editor for the text of the post, with formatting.
  - **Image**: the picture of the post.
    - **Add Picture**: chooses an image file. RetroShare then asks **Do you
      want to share and link the original image?** Answer **Yes** to share
      the original file and put a link to it in the post, so that readers can
      download it at full quality.
    - **Paste Picture** (tooltip **Paste image from clipboard**): uses the
      picture that is in the clipboard.
    - **Remove image**: removes the picture, once one is set.
    - A post is limited to 32 KB, so the picture is scaled down to fit.
      Animated pictures are not scaled; they must be under 194 KB.
  - **Link**: the **Url** field holds a web address or a RetroShare link.
    Readers open it by clicking the title of the post.
- **Post** publishes the post. Closing the window asks whether to discard it.

### Make Comment

- **Type in your comment**: the text of the comment. The dialog shows how many
  characters remain.
- **Post as:** the identity that signs the comment.
- **Post** publishes the comment.

**Reply to Comment** is the same dialog; it shows **Replying to** and the
author you answer to.

## How to

### Share a link

1. Open the board and subscribe to it if you have not done so.
2. On the **Posts** tab, click **Create Post**.
3. Type a title, click **Link** and paste the address in **Url**.
4. Choose the identity in **Post as** and click **Post**.

### Vote and comment

1. Choose your identity in the identity box of the **Posts** tab.
2. Click **Vote up** or **Vote down** on a post.
3. Click **Comments** on the post, then **Comment** in the tab that opens.

## Related settings

- [Boards](settings-boards.md): tabs.
- [People](settings-people.md): reputation rules.

## See also

- [Channels](channels.md)
- [Forums](forums.md)
- [People](people.md)
- [RetroShare links](links.md)
