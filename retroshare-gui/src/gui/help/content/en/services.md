# Optional services

RetroShare contains a few experimental services that are not part of the
standard packages. Each one adds a page to the main window, next to
[Channels](channels.md), [Forums](forums.md) and [Boards](boards.md), and is
only there when the person who built your RetroShare enabled it. If you do
not see the page, your build does not include the service.

These services are experimental: they may be incomplete, they may change, and
you can only use them with friends whose RetroShare includes them too.

Like channels and forums, they are spread from friend to friend, they are
signed with your [identities](people.md), and you receive the content of what
you subscribe to.

## Wiki

The **Wiki** page lets a group of people write and correct pages together.
Every edit of a page is kept, so the history of a page can be read and an old
version can be brought back.

It is included when RetroShare is built with the option `RS_GXSWIKIPOS`, which
is off by default.

### Layout

- On the left, the list of wiki groups, sorted in **My Groups**, **Subscribed
  Groups**, **Popular Groups** and **Other Groups**, with a **Create Group**
  button.
- On the right, the pages of the selected group in a list with the columns
  **Page Name**, **Page Id** and **Orig Id**, a **Search** field, and the text
  of the selected page.
- Above the pages, the buttons **New Page**, **Edit** and **Republish**.

### Main actions

- **Create Group**: creates a wiki group.
- **Subscribe to Group** / **Unsubscribe to Group** (context menu of a group):
  starts or stops receiving its pages. **Show Wiki Group** shows its details
  and **Edit Wiki Group** changes them.
- **New Page**: opens the **Create New Wiki Page** window. Type the **Page
  Name:** and the text, use **Preview** to see the result, and click
  **Submit**.
- **Edit** (button, or context menu of a page): opens the **Edit Wiki Page**
  window on the selected page. **History** shows the **Page Edit History**,
  where you can go back to an earlier edit with **Revert** or combine several
  edits with **Merge for Republish**.
- **Republish**: publishes the page again; offered to the people who can
  moderate the group.
- **Mark as read**, **Mark as unread**, **Show Page Id**, **Show Orig Id**:
  the other entries of the context menu of a page.

## The Wire

**The Wire** is a microblog. You create an account, called a wire, you post
short messages called pulses, and you follow the wires of other people.

It is included when RetroShare is built with the option `RS_GXSTHEWIRE`, which
is off by default.

### Layout

- At the top, the buttons **Create Account** and **Post Pulse**, and a
  **Search...** field.
- On one side, the list of wires. A box chooses which ones are listed:
  **All**, **Yourself**, **Following** or **Others**. Each wire is marked
  **Own**, **Following** or **Other** and has a **Follow** or **Unfollow**
  button.
- In the centre, the pulses. **Show Posts from** limits them to the **Last 24
  hours**, **Last 7 days**, **Last 30 days** or **All Time**.

### Main actions

- **Create Account**: opens the **Create New Wire** dialog, where you name
  your wire and can add a **Tagline:**, a **Location:** and a background
  image.
- **Follow** / **Unfollow**: subscribes to a wire or leaves it.
- **Post Pulse**: opens the pulse window. Type the text under **Whats
  happening?**, add a **URL** and up to four pictures if you want, and click
  **Post**. You need a wire of your own first.
- On a pulse, you can reply, republish it or like it. A reply can carry a
  sentiment: **Positive**, **Neutral** or **Negative**.

## Photo Albums

The **Photo Albums** page shares albums of pictures with your friends.

It is included when RetroShare is built with the option `RS_GXSPHOTOSHARE`.
It is not part of the standard builds.

### Layout

- A toolbar with **Create Album**, **View Photo**, **Slide Show** and
  **Subscribe To Album**.
- The albums, sorted in **My Albums**, **Subscribed Albums** and **Shared
  Albums**.
- The photos of the selected album.

### Main actions

- **Create Album**: opens the **Create New Album** dialog. Besides the name
  and the description, you can set a **Caption:**, the **Photographer:**,
  **Where:** and **When:**, and choose how the pictures are shared, from
  small embedded copies to the original files.
- Opening one of your albums lets you use **Add Photo**, **Edit Photo** and
  **Delete Photo**, then **Publish Photos** to send the changes.
- **Subscribe To Album**: subscribes to the selected album, so that you
  receive its photos.
- **View Photo**: opens the selected photo with its details and comments.
- **Slide Show**: shows the photos of the selected album one after the other.
  **Start** and **Stop** control it.

## See also

- [Channels](channels.md)
- [Forums](forums.md)
- [Boards](boards.md)
- [People](people.md)
