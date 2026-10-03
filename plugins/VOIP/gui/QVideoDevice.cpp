/*******************************************************************************
 * plugins/VOIP/gui/QVideoDevice.cpp                                           *
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

#include <QPainter>
#include <QImageReader>
#include <QBuffer>
#include <QTimer>
#include <QCamera>
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
#  include <QMediaDevices>
#  include <QCameraDevice>
#  include <QMediaCaptureSession>
#  include <QVideoSink>
#  include <QVideoFrameFormat>
#  include <QCameraFormat>
#else
#  include <QCameraInfo>
#endif
#if QT_VERSION >= QT_VERSION_CHECK(6, 5, 0)
#  include <QPermissions>
#  include <QApplication>
#endif
#if defined(Q_OS_MACOS)
#  include "MacCameraPermission.h"
#  include "MacCameraCapture.h"
#endif
#include <util/rsdebug.h>
#include "QVideoDevice.h"
#include "VideoProcessor.h"

// #define DEBUG_QVIDEODEVICE 1

QVideoInputDevice::QVideoInputDevice(QWidget *parent)
  :QObject(parent)
{
	_capture_device = NULL ;
	_video_processor = NULL ;
	_echo_output_device = NULL ;
	_camera_permission_granted = false ;
#if defined(Q_OS_MACOS)
    _mac_capture = NULL;
#endif
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    _capture_session = NULL;
    _video_sink = NULL;
#else
    _video_surface = NULL;
#endif
}

QVideoInputDevice::~QVideoInputDevice()
{
    stop() ;             // releases _capture_device and the capture pipeline
#if defined(Q_OS_MACOS)
    delete _mac_capture ;
    _mac_capture = NULL ;
#endif
    _video_processor = NULL ;
}

bool QVideoInputDevice::stopped() const
{
#if defined(Q_OS_MACOS)
    return _mac_capture == NULL || !_mac_capture->isRunning() ;
#else
    return _capture_device == NULL ;
#endif
}

void QVideoInputDevice::stop()
{
#if defined(Q_OS_MACOS)
    // macOS captures via AVFoundation (see start()); no Qt camera pipeline to tear
    // down. Keep the object around (cheap) so it can be restarted.
    if(_mac_capture != NULL)
        _mac_capture->stop() ;

    if(_echo_output_device != NULL)
        _echo_output_device->showFrameOff() ;
    return;
#endif

#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    _capture_device_info = QCameraDevice();

    if(_capture_device != NULL)
    {
        _capture_device->stop() ;
        delete _capture_device ;
        _capture_device = NULL ;
    }
    // The capture session wires camera -> sink; tear it down after the camera.
    delete _capture_session ;
    _capture_session = NULL ;
    delete _video_sink ;
    _video_sink = NULL ;
#else
    _capture_device_info = QCameraInfo();

	if(_capture_device != NULL)
	{
        _capture_device->stop() ;
        delete _capture_device ;     // detaches and releases the viewfinder surface
        _capture_device = NULL ;
    }
    // Delete the surface only after the camera that referenced it is gone.
    delete _video_surface ;
    _video_surface = NULL ;
#endif

    if(_echo_output_device != NULL)
        _echo_output_device->showFrameOff() ;
}
void QVideoInputDevice::getAvailableDevices(QList<QString>& device_desc)
{
    device_desc.clear();

#if defined(Q_OS_MACOS)
    MacCameraCapture::availableDevices(device_desc);
    return;
#endif

#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    const QList<QCameraDevice> dev_list = QMediaDevices::videoInputs();

    for(auto& cam:dev_list)
        device_desc.push_back(cam.description());
#else
    QList<QCameraInfo> dev_list = QCameraInfo::availableCameras();

    for(auto& cam:dev_list)
        device_desc.push_back(cam.deviceName());
#endif
}

void QVideoInputDevice::start(const QString& description)
{
	// make sure everything is re-initialised
	//
	stop() ;

#if defined(Q_OS_MACOS)
    // macOS: bypass Qt6's QPermission API entirely. qApp->requestPermission(
    // QCameraPermission) returns Denied WITHOUT ever contacting the OS -- proven with
    // the TCC log, where opening the VOIP panel triggers a kTCCServiceMicrophone
    // request (granted) but NEVER a kTCCServiceCamera request. We therefore ask
    // AVFoundation directly, which forces the real system request + prompt (same
    // mechanism the microphone already uses), then re-enter start() once granted.
    if(!_camera_permission_granted)
    {
        requestMacCameraAccess([this, description](bool granted)
        {
            if(granted)
            {
                _camera_permission_granted = true ;
                start(description) ;   // re-enter, now authorized
            }
            else
                RsDbg() << "DISTANT_VOIP: [mac] camera access not granted -- no preview" ;
        });
        return;
    }

    // Authorized: capture straight from AVFoundation. Qt's QMediaCaptureSession/
    // QVideoSink pipeline delivers zero frames on macOS (Qt 6.11), so we grab frames
    // directly; they arrive as QImage on the GUI thread and go through the very same
    // encoder + echo-display path (deliverCameraImage) as every other platform.
    if(_mac_capture == NULL)
        _mac_capture = new MacCameraCapture();

    if(_mac_capture->start(description, [this](const QImage& image){ deliverCameraImage(image); }))
        emit cameraCaptureInfo(CAMERA_IS_READY, QCamera::NoError);
    else
        emit cameraCaptureInfo(CANNOT_INITIALIZE_CAMERA, QCamera::NoError);

    return;   // macOS never falls through to the Qt camera pipeline below
#elif QT_VERSION >= QT_VERSION_CHECK(6, 5, 0)
    // Other Qt6 platforms (Linux/Windows): Qt's QPermission API works, so use it.
    // We do NOT gate on checkPermission()'s value (it can wrongly report Denied);
    // we always call requestPermission() until we get an explicit Granted.
    if(!_camera_permission_granted)
    {
        QCameraPermission cameraPermission;
        RsDbg() << "DISTANT_VOIP: camera checkPermission() status="
                << (int)qApp->checkPermission(cameraPermission) << " -- requesting anyway";

        qApp->requestPermission(cameraPermission, this,
            [this, description](const QPermission& perm)
            {
                if(perm.status() == Qt::PermissionStatus::Granted)
                {
                    _camera_permission_granted = true ;
                    start(description) ;   // re-enter, now authorized
                }
                else
                    RsDbg() << "DISTANT_VOIP: [ERROR] camera permission NOT granted after requestPermission()";
            });
        return;
    }
#endif

#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    QCameraDevice caminfo ;

    if(description.isNull())
        caminfo = QMediaDevices::defaultVideoInput();
    else
    {
        const auto cam_list = QMediaDevices::videoInputs();

        for(auto& s:cam_list)
            if(s.description() == description)
                caminfo = s;
    }

    if(caminfo.isNull())
    {
        RsDbg() << "DISTANT_VOIP: [CRITICAL] No video camera available in this system!";
        return ;
    }
    _capture_device_info = caminfo;
    RsDbg() << "DISTANT_VOIP: Initializing camera: " << caminfo.description().toStdString() << " (ID: " << caminfo.id().toStdString() << ")";

    _capture_device = new QCamera(caminfo);

    if(_capture_device->error() != QCamera::NoError)
    {
        emit cameraCaptureInfo(CANNOT_INITIALIZE_CAMERA,_capture_device->error());
        RsDbg() << "DISTANT_VOIP: [ERROR] Cannot initialise camera. Error code: " << (int)_capture_device->error();
        return;
    }

    // macOS/Qt6: QCamera can silently refuse to activate (isActive stays false,
    // errorOccurred never fires) when no capture format is selected. Choose one
    // explicitly, preferring a ~640-1280 wide format, else the first available.
    {
        const QList<QCameraFormat> formats = caminfo.videoFormats();
        RsDbg() << "DISTANT_VOIP: [Qt6] camera advertises " << formats.size() << " video format(s)";
        if(!formats.isEmpty())
        {
            QCameraFormat chosen = formats.first();
            for(const QCameraFormat& f : formats)
                if(f.resolution().width() >= 640 && f.resolution().width() <= 1280)
                {
                    chosen = f;
                    break;
                }
            _capture_device->setCameraFormat(chosen);
            RsDbg() << "DISTANT_VOIP: [Qt6] selected camera format "
                    << chosen.resolution().width() << "x" << chosen.resolution().height()
                    << " pixfmt=" << (int)chosen.pixelFormat();
        }
    }

    // Qt6 grabs frames through a QVideoSink fed by a QMediaCaptureSession
    // (QAbstractVideoSurface + setViewfinder were removed). videoFrameChanged()
    // is delivered on the GUI thread.
    _video_sink = new QVideoSink(this);
    _capture_session = new QMediaCaptureSession(this);
    _capture_session->setCamera(_capture_device);
    _capture_session->setVideoSink(_video_sink);
    // Modern pointer-to-member connect: the string-based SIGNAL/SLOT form can
    // silently mismatch on QVideoFrame. Frames are delivered on the GUI thread.
    QObject::connect(_video_sink, &QVideoSink::videoFrameChanged,
                     this, &QVideoInputDevice::handleSurfaceFrame);
    // We were blind to camera start failures: QCamera reports them asynchronously
    // via errorOccurred, NOT through error() right after construction. Without this
    // the camera can fail to activate (LED stays off, no frames) with no trace.
    QObject::connect(_capture_device, &QCamera::errorOccurred, this,
        [](QCamera::Error err, const QString& msg)
        { RsDbg() << "DISTANT_VOIP: [Qt6] QCamera errorOccurred=" << (int)err
                  << " : " << msg.toStdString(); });
    QObject::connect(_capture_device, &QCamera::activeChanged, this,
        [](bool active)
        { RsDbg() << "DISTANT_VOIP: [Qt6] QCamera activeChanged=" << active; });
    RsDbg() << "DISTANT_VOIP: Capture session wired camera -> video sink.";
#else
    QCameraInfo caminfo ;

    if(description.isNull())
        caminfo = QCameraInfo::defaultCamera();
    else
    {
        auto cam_list = QCameraInfo::availableCameras();

        for(auto& s:cam_list)
            if(s.deviceName() == description)
                caminfo = s;
    }

    if(caminfo.isNull())
    {
        RsDbg() << "DISTANT_VOIP: [CRITICAL] No video camera available in this system!";
        return ;
    }
    _capture_device_info = caminfo;
    RsDbg() << "DISTANT_VOIP: Initializing camera: " << caminfo.description().toStdString() << " (ID: " << caminfo.deviceName().toStdString() << ")";

    _capture_device = new QCamera(caminfo);

    if(_capture_device->error() != QCamera::NoError)
    {
        emit cameraCaptureInfo(CANNOT_INITIALIZE_CAMERA,_capture_device->error());
        RsDbg() << "DISTANT_VOIP: [ERROR] Cannot initialise camera. Error code: " << (int)_capture_device->error();
        return;
    }

    _capture_device->setCaptureMode(QCamera::CaptureVideo);

    // Grab frames through a viewfinder surface. This works on every backend,
    // unlike QVideoProbe which attaches but never delivers frames on macOS
    // (AVFoundation) -> camera LED on, but no image.
    _video_surface = new RsCameraVideoSurface(this);
    QObject::connect(_video_surface,SIGNAL(frameAvailable(QVideoFrame)),this,SLOT(handleSurfaceFrame(QVideoFrame)));
    _capture_device->setViewfinder(_video_surface);
    RsDbg() << "DISTANT_VOIP: Viewfinder surface attached to camera.";
#endif

    QObject::connect(this,SIGNAL(cameraCaptureInfo(CameraStatus,QCamera::Error)),this,SLOT(errorHandling(CameraStatus,QCamera::Error)));

    if(_capture_device->error() == QCamera::NoError)
    {
        RsDbg() << "DISTANT_VOIP: Camera object created.";
        emit cameraCaptureInfo(CAMERA_IS_READY,QCamera::NoError);
    }

    RsDbg() << "DISTANT_VOIP: Finalizing camera start().";
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    // Qt6/macOS: starting the camera inline, in the same block that just wired the
    // QMediaCaptureSession, can leave it inactive (isActive stays false, no error,
    // activeChanged never fires). Defer start() to the next event-loop iteration so
    // the session is fully attached first, then re-check activation once more later
    // to tell "starts late" apart from "never starts".
    {
        QCamera* cam = _capture_device;
        QTimer::singleShot(0, this, [this, cam]()
        {
            if(cam != _capture_device) return;   // stopped/replaced meanwhile
            _capture_device->setActive(true);
            RsDbg() << "DISTANT_VOIP: [Qt6] deferred setActive(true): isActive="
                    << _capture_device->isActive() << " error=" << (int)_capture_device->error();
            QTimer::singleShot(1500, this, [this, cam]()
            {
                if(cam != _capture_device) return;
                RsDbg() << "DISTANT_VOIP: [Qt6] +1.5s isActive=" << _capture_device->isActive();
            });
        });
    }
#else
    RsDbg() << "DISTANT_VOIP: LED should turn on now.";
    _capture_device->start();
#endif
}

void QVideoInputDevice::handleSurfaceFrame(const QVideoFrame& f)
{
    static int p_id = 0;
    grabFrame(p_id++, f);
}

void QVideoInputDevice::errorHandling(CameraStatus status,QCamera::Error error)
{
#ifdef DEBUG_QVIDEODEVICE
    std::cerr << "Received msg from camera capture: status=" << (int)status << " error=" << (int)error << std::endl;
#else
    Q_UNUSED(error);
#endif
    if(status == CANNOT_INITIALIZE_CAMERA)
    {
        std::cerr << "Cannot initialize camera. Make sure a QtMultimedia camera backend/plugin is installed, as this is a common cause for the camera not being found." << std::endl;
    }
}

void QVideoInputDevice::deliverCameraImage(const QImage& imageIn)
{
    // macOS AVFoundation path: frames already arrive as QImage on the GUI thread.
    // Mirror grabFrame()'s policy: keep VOIP frames modest (<= 640 wide), feed the
    // encoder, and echo a mirrored self-view.
    if(imageIn.isNull())
        return;

    QImage image = (imageIn.width() > 640) ? imageIn.scaled(640, 480, Qt::KeepAspectRatio) : imageIn;

    static int mac_frame_count = 0;
    if(mac_frame_count++ % 20 == 0)
        RsDbg() << "DISTANT_VOIP: [mac] frame delivered. Size=" << image.width() << "x" << image.height();

    if(_video_processor != NULL)
    {
        _video_processor->processImage(image) ;
        emit networkPacketReady() ;
    }
    if(_echo_output_device != NULL)
        // Self-view mirrored (like every mainstream video app); the frame sent to
        // the peer (processImage above) stays un-mirrored.
        _echo_output_device->showFrame(image.mirrored(true, false)) ;
}

void QVideoInputDevice::grabFrame(int id,QVideoFrame frame)
{
    if(frame.size().isEmpty())
    {
        RsDbg() << "DISTANT_VOIP: [WARNING] Received an empty frame from camera!";
        return;
    }

    static int frame_count = 0;
    if (frame_count++ % 20 == 0) {
        RsDbg() << "DISTANT_VOIP: Frame received. ID=" << id << " Size=" << frame.width() << "x" << frame.height() << " Format=" << (int)frame.pixelFormat();
    }

    QImage image;

#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    frame.map(QVideoFrame::ReadOnly);

    // Check if it's already a JPEG (common on some Windows cams)
    if (frame.pixelFormat() == QVideoFrameFormat::Format_Jpeg) {
        QByteArray data((const char *)frame.bits(0), frame.mappedBytes(0));
        QBuffer buffer;
        buffer.setData(data);
        buffer.open(QIODevice::ReadOnly);
        QImageReader reader(&buffer, "JPG");
        reader.setScaledSize(QSize(640,480));
        image = reader.read();
    } else {
        // Qt6's QVideoFrame::toImage() handles YUV/RGB conversion robustly.
        image = frame.toImage();

        if (!image.isNull() && (image.width() > 640))
            image = image.scaled(640, 480, Qt::KeepAspectRatio);
    }
#else
    frame.map(QAbstractVideoBuffer::ReadOnly);

    // Check if it's already a JPEG (common on some Windows cams)
    if (frame.pixelFormat() == QVideoFrame::Format_Jpeg) {
        QByteArray data((const char *)frame.bits(), frame.mappedBytes());
        QBuffer buffer;
        buffer.setData(data);
        buffer.open(QIODevice::ReadOnly);
        QImageReader reader(&buffer, "JPG");
        reader.setScaledSize(QSize(640,480));
        image = reader.read();
    } else {
        // Use Qt's native conversion which handles YUV, RGB, etc.
        // If your Qt version is 5.15+, frame.image() is available.
        // Otherwise, we create a QImage from the raw bits.
        image = frame.image();

        if (image.isNull()) {
            // Manual fallback for common formats if .image() fails
            image = QImage(frame.bits(), frame.width(), frame.height(), frame.bytesPerLine(), QVideoFrame::imageFormatFromPixelFormat(frame.pixelFormat()));
        }

        if (!image.isNull() && (image.width() > 640)) {
            image = image.scaled(640, 480, Qt::KeepAspectRatio);
        }
    }
#endif

    if (image.isNull()) {
        RsDbg() << "DISTANT_VOIP: [ERROR] Failed to convert QVideoFrame to QImage. Format=" << (int)frame.pixelFormat();
        frame.unmap();
        return;
    }

    frame.unmap();

#ifdef DEBUG_QVIDEODEVICE
    std::cerr << "Frame " << id << ". Pixel format: " << frame.pixelFormat() << ". Size: " << image.size().width() << " x " << image.size().height() << std::endl; // if(frame.pixelFormat() != QVideoFrame::Format_Jpeg)
#else
    Q_UNUSED(id);
#endif

    if(_video_processor != NULL)
    {
        _video_processor->processImage(image) ;

        emit networkPacketReady() ;
    }
    if(_echo_output_device != NULL)
        // Self-view is mirrored (horizontal flip), like every mainstream video
        // app, so it feels natural. The frame sent to the peer (processImage
        // above) stays un-mirrored, so the remote sees us the right way round.
        _echo_output_device->showFrame(image.mirrored(true, false)) ;
}

bool QVideoInputDevice::getNextEncodedPacket(RsVOIPDataChunk& chunk)
{
    if(_video_processor)
        return _video_processor->nextEncodedPacket(chunk) ;
    else
        return false ;
}

uint32_t QVideoInputDevice::currentBandwidth() const
{
    if(stopped())
        return 0;
    else
        return _video_processor->currentBandwidthOut() ;
}

QVideoOutputDevice::QVideoOutputDevice(QWidget *parent)
    : QLabel(parent)
{
	showFrameOff() ;
}

void QVideoOutputDevice::showFrameOff()
{
	setPixmap(QPixmap(":/images/video-icon-big.png").scaled(QSize(height()*4/3,height()),Qt::KeepAspectRatio,Qt::SmoothTransformation)) ;
	setAlignment(Qt::AlignCenter);
}

void QVideoOutputDevice::showFrame(const QImage& img)
{
#ifdef DEBUG_QVIDEODEVICE
    std::cerr << "img.size = " << img.width() << " x " << img.height() << std::endl;
#endif
    setPixmap(QPixmap::fromImage(img).scaled( QSize(height()*4/3,height()),Qt::KeepAspectRatio,Qt::SmoothTransformation)) ;
}
