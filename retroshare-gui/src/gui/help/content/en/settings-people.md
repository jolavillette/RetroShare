# People

The People page of the [Preferences](settings.md) sets how the reputation of
identities is computed and how long banned identities are kept.

An identity has an opinion that you give it yourself (negative, neutral or
positive) and an opinion that comes from the votes of your friend nodes. A
banned identity is one with a negative reputation: its posts and messages are
not shown or forwarded. See [People](people.md).

## Buttons

All the options are in the **Reputation** group.

- **use "positive" as the default opinion for contacts (instead of neutral)**:
  an identity that is in your contacts and to which you gave no opinion gets
  a positive opinion automatically. This raises the reputation of the
  identities you actually use. On by default.
- **Automatically add identities owned by friend nodes to my contacts**: the
  identities signed by the profile of one of your friends are added to your
  contacts without your doing anything. On by default.
- **Difference in votes (+/-) to rate an ID positively:** the number of
  positive votes of your friends that an identity needs above its negative
  votes to go from "Neutral" to "Positive (according to your friends)". From
  1 to 100, 1 by default. A higher number makes it harder.
- **Difference in votes (+/-) to rate an ID negatively:** the number of
  negative votes of your friends that an identity needs above its positive
  votes to go from "Neutral" to "Negative (according to your friends)". From
  1 to 100, 1 by default.
- **Delete banned identities after (0 means never):** in days. A banned
  identity that has not been used for this time is deleted from your node.
  With 0, banned identities are never deleted for that reason.
- **Reset reputation of banned identities after (0 means never):** in days.
  After you have banned an identity, RetroShare remembers the ban for this
  time, so that the identity does not come back when a forum or a channel
  uses it. After the delay the ban is cleared, and the identity is
  downloaded again as not banned if something uses it. With 0, the ban is
  never cleared.

## See also

- [People](people.md)
- [Glossary](glossary.md)
