/*
-----------------------------------------------------------------------------
Filename:    BaseApplication.h
-----------------------------------------------------------------------------

This source file is part of the
   ___                 __    __ _ _    _
  /___\__ _ _ __ ___  / / /\ \ (_) | _(_)
 //  // _` | '__/ _ \ \ \/  \/ / | |/ / |
/ \_// (_| | | |  __/  \  /\  /| |   <| |
\___/ \__, |_|  \___|   \/  \/ |_|_|\_\_|
      |___/
      Tutorial Framework
      http://www.ogre3d.org/tikiwiki/
-----------------------------------------------------------------------------
*/
#ifndef __BaseApplication_h_
#define __BaseApplication_h_

#include "utils/Utilities.h"
DISABLE_WARNINGS

#include <OgreCamera.h>
#include <OgreConfigFile.h>
#include <OgreEntity.h>
#include <OgreLogManager.h>
#include <OgreRenderWindow.h>
#include <OgreRoot.h>
#include <OgreSceneManager.h>
#include <OgreViewport.h>

#include <OISEvents.h>
#include <OISInputManager.h>
#include <OISKeyboard.h>
#include <OISMouse.h>

#include <SdkCameraMan.h>
#include <SdkTrays.h>

#ifdef OGRE_STATIC_LIB
#define OGRE_STATIC_GL
#define OGRE_STATIC_ParticleFX
#define OGRE_STATIC_BSPSceneManager
//#define OGRE_STATIC_CgProgramManager
#define OGRE_STATIC_OctreeZone
#define OGRE_STATIC_OctreeSceneManager
#include "OgreStaticPluginLoader.h"
#endif

ENABLE_WARNINGS

class BaseApplication : public Ogre::FrameListener,
                        public Ogre::WindowEventListener,
                        public OIS::KeyListener,
                        public OIS::MouseListener,
                        OgreBites::SdkTrayListener
{
public:
    BaseApplication( void );
    virtual ~BaseApplication( void );

    virtual void go( void );

protected:
    virtual bool setup();
    virtual bool configure( void );
    virtual void chooseSceneManager( void );
    virtual void createCamera( void );
    virtual void createFrameListener( void );
    virtual void createScene( void ) = 0; // Override me!
    virtual void destroyScene( void );
    virtual void createViewports( void );
    virtual void setupResources( void );
    virtual void createResourceListener( void );
    virtual void loadResources( void );

    // Ogre::FrameListener
    virtual bool frameRenderingQueued( const Ogre::FrameEvent &evt );

    // OIS::KeyListener
    virtual bool keyPressed( const OIS::KeyEvent &arg );
    virtual bool keyReleased( const OIS::KeyEvent &arg );
    // OIS::MouseListener
    virtual bool mouseMoved( const OIS::MouseEvent &arg );
    virtual bool mousePressed( const OIS::MouseEvent &arg, OIS::MouseButtonID id );
    virtual bool mouseReleased( const OIS::MouseEvent &arg, OIS::MouseButtonID id );

    // Ogre::WindowEventListener
    // Adjust mouse clipping area
    virtual void windowResized( Ogre::RenderWindow *rw );
    // Unattach OIS before window shutdown (very important under Linux)
    virtual void windowClosed( Ogre::RenderWindow *rw );

    Ogre::Root *mRoot;
    Ogre::Camera *mCamera;
    Ogre::SceneManager *mSceneMgr;
    Ogre::RenderWindow *mWindow;
    Ogre::String mResourcesCfg;
    Ogre::String mPluginsCfg;
    Ogre::OverlaySystem *mOverlaySystem;

    // OgreBites
    OgreBites::SdkTrayManager *mTrayMgr;
    OgreBites::SdkCameraMan *mCameraMan;   // basic camera controller
    OgreBites::ParamsPanel *mDetailsPanel; // sample details panel
    bool mCursorWasVisible;                // was cursor visible before dialog appeared
    bool mShutDown;

    // OIS Input devices
    OIS::InputManager *mInputManager;
    OIS::Mouse *mMouse;
    OIS::Keyboard *mKeyboard;

#ifdef OGRE_STATIC_LIB
    Ogre::StaticPluginLoader m_StaticPluginLoader;
#endif
};

#endif // #ifndef __BaseApplication_h_
