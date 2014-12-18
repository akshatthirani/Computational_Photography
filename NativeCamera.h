// TAGRELEASE: PUBLIC
#ifndef NativeCamera_H
#define NativeCamera_H

#include "/Users/Akshat/NVPACK/Samples/GameWorks_Samples/extensions/include/NvAppBase/NvSampleApp.h"
#include "/Users/Akshat/NVPACK/Samples/GameWorks_Samples/extensions/include/NvGLUtils/NvImage.h"

#include "/Users/Akshat/Documents/Eclipse Workspace hw2/HW2Incomplete/native_camera2/include/native_camera2/native_camera2.h"

#include <memory>
#include <sstream>
#include <string>

#include <iostream>
#include <chrono>
#include <thread>

class NvStopWatch;
class NvFramerateCounter;

class NativeCamera: public NvSampleApp {
public:
	NativeCamera(NvPlatformContext* platform);
	~NativeCamera();

	void initUI(void);
	void initRendering(void);
	void focusChanged(bool focused) override;
	void draw(void);
	void reshape(int32_t width, int32_t height);

	void configurationCallback(NvEGLConfiguration& config);

protected:
	void drawStreamImage(const nv::camera2::CameraStream *stream,
			GLuint *textures);
	void setupStreamTextures(const nv::camera2::CameraStream *stream,
			GLuint *textures);
	void uploadImage(nv::camera2::CameraBuffer &img, GLuint *textures);

	void startCamera();
	void stopCamera();
	std::unique_ptr<nv::camera2::CameraStream> setupCameraStream(
			nv::camera2::PixelFormat format);
	void capture();

	std::unique_ptr<NvFramerateCounter> m_framerate;
	std::unique_ptr<NvGLSLProgram> m_progYUV;
	std::unique_ptr<NvGLSLProgram> m_progRAW;
	GLuint m_imgTextures[4];
	float m_viewAspectRatio;
	float m_imgAspectRatio;

	// Camera objects
	std::unique_ptr<nv::camera2::CameraManager> m_cameraManager;
	std::unique_ptr<nv::camera2::CameraDevice> m_cameraDevice;
	std::unique_ptr<nv::camera2::CameraStream> m_streamYUV;
	std::unique_ptr<nv::camera2::CameraStream> m_streamRAW;
	std::unique_ptr<nv::camera2::CameraStream> m_streamJPG;

	nv::camera2::CaptureRequest req;

	nv::camera2::StaticProperties m_staticProperties;
	nv::camera2::CaptureRequest m_request;
	nv::camera2::Size m_imgSize;

	// TweakBar
	bool m_useRawCapture = false;
	bool m_AELock = false;
	bool m_takeShot = false;
	float m_AECompensation = 0.0f;
	int test = 0;
	int numShots = 0;
	std::ostringstream stm;
	std::string fName;
	uint32_t mSensitivity;
	uint32_t mExposureTime;
	uint32_t mFocusDistance;


};

#endif

