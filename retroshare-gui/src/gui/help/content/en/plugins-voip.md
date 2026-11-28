# VOIP

The VOIP plugin adds audio and video calls to RetroShare. A call takes place
in the private chat window that you open with a friend node, and goes
directly over the encrypted connection between your two nodes.

Conditions:

- VOIP is a [plugin](plugins.md): it must be present and enabled in the
  [Plugins preferences](settings-plugins.md), on your side and on your
  friend's side.
- You can only call a friend node that is connected. The call buttons are
  added to the private chat with a friend node; they do not exist in distant
  chats (with an identity that is not your friend) nor in chat rooms.
- You need a microphone and loudspeakers or a headset, and a camera for
  video. RetroShare uses the default sound input, sound output and camera of
  your system.

The plugin is experimental.

## Layout

In the title bar of a private chat window, the plugin adds six buttons. They
carry icons only; hover them to read their current tooltip.

- **Mute** (loudspeaker): only available during an audio call.
- **Start Call** (microphone): starts and stops the audio.
- **Start Video Call** (camera): starts and stops the video.
- **Hangup Call** (red telephone): shown only during a call.
- **Hide Chat Text** (bubble): only available during a video call.
- **Fullscreen mode**: only available during a video call.

Two small progress indicators appear next to the buttons while a call is
ringing, one for audio and one for video.

During a video call, the picture of your friend is shown above the chat text,
with the picture of your own camera in small.

## Buttons

- **Start Call**: calls your friend for an audio conversation. While it
  rings, the audio indicator turns and its tooltip reads **Waiting for your
  friend to respond to your audio call.** When your friend accepts, the chat
  shows that they accepted your audio call, the sound starts, and the tooltip
  of the button becomes **Hold Call**. Click the button again to stop your
  audio; the tooltip becomes **Resume Call**. Clicking the button while your
  own call is still ringing cancels it.
- **Start Video Call**: calls your friend for a video conversation, in the
  same way. Once accepted, the chat shows **You're now sending video...**,
  the pictures appear, and the tooltip of the button becomes **Shut camera
  off**. Click again to stop the video; the tooltip becomes **Activate
  camera**. Audio and video are two separate calls: start both if you want to
  see and hear each other.
- **Hangup Call**: ends the audio and the video at once.
- **Hide Chat Text** / **Show Chat Text**: hides the text of the conversation
  to leave the room to the video, and shows it again.
- **Fullscreen mode** / **Return to normal view.**: shows the video on the
  whole screen, with a copy of the buttons at the bottom of the screen, and
  comes back to the chat window.
- **Mute**: the button is pressed in automatically when the audio starts and
  its tooltip then alternates between **Mute yourself** and **Unmute
  yourself**. In this version it does not change the sound. To stop sending
  your voice, use **Hold Call**.

Every step of a call is written in the chat as a **VoIP Status** line:
**Outgoing Call is started...**, **Outgoing Audio Call stopped.**, **Video
call stopped**, or that your friend hung up.

## Receiving a call

When a friend calls you:

- a ringing sound is played once per second;
- a toaster appears with the name of your friend followed by "calling" and
  two buttons: **Answer** (or **Answer with video** for a video call) accepts
  the call and opens the chat window, **Decline** refuses it;
- in the private chat window with this friend, a **VoIP Status** line says
  that your friend invites you to an audio or a video conversation, with two
  buttons in the text: **Accept Audio Call** and **Decline Audio Call**, or
  **Accept Video Call** and **Decline Video Call**;
- the microphone or camera button of the title bar blinks and its tooltip
  reads **Answer**: clicking it accepts the call as well.

When you accept, your microphone (or your camera) is turned on and the call
starts. When you decline, your friend's call stops.

## The VOIP preferences page

The plugin adds a **VOIP** page to the preferences. While this page is
displayed, your microphone and your camera are on, so that you can test them.
Nothing is sent to anyone. The audio changes apply immediately.

### Video

- The picture on the left shows what your camera captures.
- **Input device:** chooses the camera used for this preview, or **[No
  video]**. The list is refreshed when you plug a camera in. During a call
  RetroShare uses the default camera of your system.
- **Available bandwidth:** (in **KB/s**): simulates the bandwidth available
  for the video, so that you can judge the picture quality at that rate. It
  does not limit real calls.
- **preview**: when checked, the picture shows the video as your friend
  would receive it, after compression at the simulated rate. When unchecked,
  it shows the raw picture of the camera.
- The graph **Required bandwidth** draws the rate that the compressed video
  currently needs.

During a call, the video rate adapts by itself to what your friend's node
reports it can receive.

### Audio

- The bar at the top shows the level of your microphone. With **Voice
  Activity** it is split into three zones: silence, uncertain, speech.
- **Noise Suppression**: how strongly steady background noise (fan, hum) is
  removed. The value is shown in dB, or **Off**.
- **Amplification**: the maximum amplification applied to your microphone.
  RetroShare brings your voice to a normal volume before compressing it; this
  limits how much it may amplify.
- **Echo Cancellation Processing**: removes from your microphone the sound
  that comes out of your own loudspeakers, so that your friend does not hear
  themselves. Leave it on unless you use a headset.
- **Transmit:** when your voice is sent:
    - **Continuous**: all the time during the call;
    - **Voice Activity**: only when RetroShare detects that you speak. Three
      more settings appear, see below;
    - **Push To Talk**: do not choose it. This version has no key to press to
      talk, so nothing is transmitted in this mode.
- **Audio Wizard**: opens the wizard described below.

With **Voice Activity**:

- **Voice Hold**: how long the transmission goes on after you stop speaking.
  Raise it if your words are cut.
- **Silence Below**: microphone levels under this value count as silence.
- **Speech Above**: levels above this value count as speech. Between the two
  values, transmission goes on if you were already speaking.

## Dialogs

### Audio Tuning Wizard

The wizard helps you set the microphone in four steps. During the wizard you
hear your own microphone played back.

1. **Introduction**.
2. **Volume tuning**: open the sound settings of your system, select the
   microphone as input with the maximum recording volume, then talk normally
   and move the slider so that the bar reaches the green zone when you talk
   and does not enter the orange zone. **Stop looping echo for this wizard**
   stops the playback of your own voice. **Use high contrast graphics**
   changes the colours of the bars for visually impaired users.
3. **Voice Activity Detection**: choose **Voice Detection** and move the
   slider so that your first words fall in the green zone, your speech stays
   in the yellow zone and silence stays in the red zone; or choose
   **Continuous transmission**. Leave **Push To Talk:** aside.
4. **Finished**.

## Notifications

- **Toasters**: in the [Notify preferences](settings-notify.md), tab
  **Toasters**, the group **VOIP** has two boxes, **Audio Call** and **Video
  Call**, to show or not the toaster of an incoming call.
- **Sounds**: in the [Sound preferences](settings-sound.md), the group
  **VOIP** has four events: **Incoming audio call**, **Incoming video call**,
  **Outgoing audio call** and **Outgoing video call**.

## How to

### Check your microphone before a first call

1. Open [Preferences](settings.md) and select **VOIP**.
2. Speak: the bar must move. If it does not, check the default input device
   in the sound settings of your system.
3. Click **Audio Wizard** and follow the steps.
4. Keep **Transmit:** on **Voice Activity** or **Continuous**.

### Call a friend

1. In the [Network](network.md) page, open the private chat with a friend
   node that is connected.
2. Click the microphone button (**Start Call**) and wait for your friend to
   answer.
3. Click the camera button (**Start Video Call**) if you also want video.
4. Click the red telephone (**Hangup Call**) to end the call.

## Related settings

- [Plugins](settings-plugins.md): enabling the VOIP plugin.
- [Notify](settings-notify.md): toasters of incoming calls.
- [Sound](settings-sound.md): ring tones.
- [Chats](settings-chat.md): the private chat window.

## See also

- [Plugins](plugins.md)
- [Chats](chats.md)
- [Network](network.md)
