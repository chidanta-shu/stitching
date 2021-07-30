// dllmain.cpp : 定义 DLL 应用程序的入口点。
#include "stdafx.h"
//#ifdef LOG_OUTPUT
//#include <glog/logging.h>
//#endif

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
	char    szBuff[MAX_PATH] = { 0 };
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
//#ifdef LOG_OUTPUT
//		GetModuleFileNameA(hModule, szBuff, MAX_PATH);
//
//		FLAGS_log_dir = "./log";
//		// 初始化GLog库 		
//		google::InitGoogleLogging(szBuff); 
//		google::SetLogDestination(google::GLOG_INFO, "./log/INFO_");
//		google::SetLogDestination(google::GLOG_WARNING, "./log/WARNING_");
//		google::SetLogDestination(google::GLOG_ERROR, "./log/ERROR_");
//
//		google::SetStderrLogging(google::GLOG_FATAL);
//		google::SetLogFilenameExtension("log_");
//		FLAGS_colorlogtostderr = true;  // Set log color
//		FLAGS_logbufsecs = 0;  // Set log output speed(s)
//		FLAGS_max_log_size = 1024;  // Set max log file size
//		FLAGS_stop_logging_if_full_disk = true;  // If disk is full
//		LOG(WARNING) << "start!";
//#endif
		break;
	case DLL_THREAD_ATTACH:
		break;
	case DLL_THREAD_DETACH:
		break;
	case DLL_PROCESS_DETACH:
//#ifdef LOG_OUTPUT
//		google::ShutdownGoogleLogging(); 
//#endif
		break;
	}
	return TRUE;
}

