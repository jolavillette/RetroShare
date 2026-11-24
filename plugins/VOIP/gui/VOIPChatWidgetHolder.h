/*******************************************************************************
 * plugins/VOIP/gui/VOIPChatWidgetHolder.h                                     *
 *                                                                             *
 * Copyright (C) 2012 by Retroshare Team <retroshare.project@gmail.com>        *
 *                                                                             *
 * This program is free software: you can redistribute it and/or modify        *
 * it under the terms of the GNU Affero General Public License as              *
 * published by the Free Software Foundation, either version 3 of the          *
 * License, or (at your option) any later version.                             *
 *                                                                             *
 * This program is distributed in the hope that it will be useful,             *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of              *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                *
 * GNU Affero General Public License for more details.                         *
 *                                                                             *
 * You should have received a copy of the GNU Affero General Public License    *
 * along with this program. If not, see <https://www.gnu.org/licenses/>.       *
 *                                                                             *
 *******************************************************************************/

#pragma once
//Qt
#include <QObject>
#include <QGraphicsEffect>
#include <QProgressBar>
#include <QTimer>
#include <QElapsedTimer>
//VOIP
#include "gui/qtmultimedia_compat.h"
#include "gui/VOIPNotify.h"
#include <gui/SpeexProcessor.h>
#include "services/rsVOIPItems.h"
//retroshare-gui
#include <gui/chat/ChatWidget.h>
#include <gui/common/RsButtonOnText.h>

class QToolButton;
class QVideoInputDevice ;
class QVideoOutputDevice ;
class VideoProcessor ;

#define VOIP_SOUND_INCOMING_AUDIO_CALL "VOIP_incoming_audio_call"
#define VOIP_SOUND_INCOMING_VIDEO_CALL "VOIP_incoming_video_call"
#define VOIP_SOUND_OUTGOING_AUDIO_CALL "VOIP_outgoing_audio_call"
#define VOIP_SOUND_OUTGOING_VIDEO_CALL "VOIP_outgoing_video_call"

class VOIPChatWidgetHolder : public QObject, public ChatWidgetHolder
{
	Q_OBJECT

public:
	VOIPChatWidgetHolder(ChatWidget *chatWidget, VOIPNotify *notify);
	virtual ~VOIPChatWidgetHolder();

	virtual void updateStatus(RsStatusValue status);

	void addAudioData(const RsPeerId &peer_id, QByteArray* array) ;
	void addVideoData(const RsPeerId &peer_id, QByteArray* array) ;
	void setAcceptedBandwidth(uint32_t bytes_per_sec) ;

	void ReceivedInvitation(const RsPeerId &peer_id, int flags) ;
	void ReceivedVoipHangUp(const RsPeerId &peer_id, int flags) ;
	void ReceivedVoipAccept(const RsPeerId &peer_id, int flags) ;

public slots:
	void sendAudioData();
	void sendVideoData();
	void startAudioCapture();
	void startVideoCapture();
	void hangupCallAudio() ;
	void hangupCallVideo() ;

private slots:
	void showCameraProblem(const QString& msg);
	void audioInputCheckTimeOut();
	void toggleMicMute();
	void toggleMicMuteFS();
	void toggleCamera();
	void toggleCameraFS();
	// Decline an incoming call (red button while it rings).
	void declineAudioCall();
	void declineVideoCall();
	// Red button: declines a call ringing in, else hangs up.
	void hangupOrDecline();
	// Answers a video call with our camera off.
	void answerWithoutCamera();
	// Call buttons: start a call, or answer one that rings in.
	void audioCallButtonClicked();
	void videoCallButtonClicked();
	void toggleAudioListen();
	void toggleAudioListenFS();
	void toggleAudioCapture();
	void toggleAudioCaptureFS();
	void toggleVideoCapture();
	void toggleVideoCaptureFS();
	void toggleHideChatText();
	void toggleFullScreen();
	void toggleFullScreenFS();
	void hangupCall() ;
	void botMouseEnterTake();
	void botMouseLeaveTake();
	void botMouseEnterDecline();
	void botMouseLeaveDecline();
	void timerAudioRingTimeOut();
	void timerVideoRingTimeOut();
	void checkTunnelAndAutoCall();

private:
	RsPeerId getEffectivePeerId(bool allowPopup = false);
	// Id the VOIP notifications know the correspondent by, without any tunnel
	// lookup: the friend in a direct chat, the tunnel (virtual peer) in a
	// distant one.
	RsPeerId notifyPeerId();
	void deleteButtonMap(int flags = RS_VOIP_FLAGS_AUDIO_DATA | RS_VOIP_FLAGS_VIDEO_DATA);
	void addNewVideoButtonMap(const RsPeerId &peer_id);
	void addNewAudioButtonMap(const RsPeerId &peer_id);

	QTimer* mTunnelWaitTimer;
	RsChats* mChats; // already exists? no, we use global, that's fine.
	void replaceFullscreenWidget();
	void showNormalView();

protected:
	bool eventFilter(QObject *obj, QEvent *event);

	// Audio input/output
	RsAudioInput* inputAudioDevice;
	RsAudioOutput* outputAudioDevice;

	QtSpeex::SpeexInputProcessor* inputAudioProcessor;
	QtSpeex::SpeexOutputProcessor* outputAudioProcessor;

	// Video input/output
	QVideoOutputDevice *outputVideoDevice;
	QVideoOutputDevice *echoVideoDevice;
	QVideoInputDevice *inputVideoDevice;

	//For FullScreen Mode
	QFrame *fullScreenFrame;
	QVideoOutputDevice *outputVideoDeviceFS;
	QVideoOutputDevice *echoVideoDeviceFS;
	Qt::WindowFlags outputVideoDeviceFlags;

	QWidget *videoWidget ;	// pointer to call show/hide

	VideoProcessor *videoProcessor;

	// Additional buttons to the chat bar
	QToolButton *audioListenToggleButton ;
	QToolButton *audioListenToggleButtonFS ;
	QToolButton *audioCaptureToggleButton ;
	QToolButton *audioCaptureToggleButtonFS ;
	QToolButton *videoCaptureToggleButton ;
	QToolButton *videoCaptureToggleButtonFS ;
	QToolButton *hideChatTextToggleButton ;
	QToolButton *fullscreenToggleButton ;
	QToolButton *fullscreenToggleButtonFS ;
	QToolButton *hangupButton ;
	QToolButton *hangupButtonFS ;
	QFrame *toolBarFS;

	typedef QMap<QString, QPair<RSButtonOnText*, RSButtonOnText*> > button_map;
	button_map buttonMapTakeCall;

	//Waiting for peer accept
	QProgressBar *pbAudioRing;
	QProgressBar *pbVideoRing;
	QTimer *timerAudioRing;
	QTimer *timerVideoRing;
	int sendAudioRingTime; //(-2 connected, -1 reseted, >=0 in progress)
	int sendVideoRingTime; //(-2 connected, -1 reseted, >=0 in progress)
	int recAudioRingTime; //(-2 connected, -1 reseted, >=0 in progress)
	int recVideoRingTime; //(-2 connected, -1 reseted, >=0 in progress)

	//TODO, remove this when soundManager can manage multi events.
	int lastTimePlayOccurs;

	VOIPNotify *mVOIPNotify;

private:
	// Device problems go to the chat as VoIP status lines, once per call start:
	// a device that cannot be opened, or opens and delivers nothing.
	void addVoipStatusMessage(const QString& text);
	// Name of the correspondent for the call messages.
	QString peerName(const RsPeerId& peer_id = RsPeerId());
	void checkAudioDevicesAtCallStart();
	void showAudioInputProblem(const QString& text);

	QTimer *timerAudioInputCheck;
	qint64 audioInputBytesAtStart;
	bool audioInputProblemShown;
	bool audioOutputProblemShown;

	// At rest the window shows the two call buttons; while a call rings or
	// runs, a call bar replaces them (microphone, camera, full screen,
	// hang-up). Shows the buttons of the current state.
	void updateCallButtons();
	// A video call carries sound: its microphone and speakers start and stop
	// with it, without an audio ring of their own.
	void videoCallStarted();
	void videoCallEnded();

	// Call bar: mute the microphone, turn the camera off, both local, the call
	// goes on.
	QToolButton *micMuteButton;
	QToolButton *micMuteButtonFS;
	QToolButton *cameraToggleButton;
	QToolButton *cameraToggleButtonFS;
	bool audioStartedWithVideo;

	// Shown while a video call rings in, next to the green video button.
	QToolButton *answerNoCameraButton;
	bool startWithCameraOff;

	// A call rings at most this long, in or out (see the ring timers).
	QElapsedTimer ringClock;
};

