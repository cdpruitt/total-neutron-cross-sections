/******************************************************************************
*
* CAEN SpA - Front End Division
* Via Vetraia, 11 - 55049 - Viareggio ITALY
* +390594388398 - www.caen.it
*
***************************************************************************//**
* \note TERMS OF USE:
* This program is free software; you can redistribute it and/or modify it under
* the terms of the GNU General Public License as published by the Free Software
* Foundation. This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. The user relies on the
* software, documentation and results solely at his own risk.
******************************************************************************/
/**
* @author Ron Fox <ron@caentech.com>
* @brief implementation of named pipe functions.
*
*/

#ifdef WIN32
#include "pipefunctions.h"
#include <stdio.h>
#include <windows.h>

/**
* CreateServerPipe
*    Creates a named pipe for the server process. The pipe is made as
*    a unidirectional conduit from the server to the single client
*    allowed.
* 
* @param basename - the base name of the pipe.  The actual pipe name will
*                   be \\.\pipe\<basename>
* @return HANDLE - handle to the pipe or INVALID_HANDLE_VALUE if there is an error.
*/
HANDLE CreateServerPipe(const char* basename)
{
	wchar_t pipename[MAX_PATH];
	HANDLE pipeHandle;

	swprintf(pipename, MAX_PATH, L"\\\\.\\pipe\\%s", basename);

	pipeHandle = CreateNamedPipe(
		pipename, PIPE_ACCESS_OUTBOUND | FILE_FLAG_FIRST_PIPE_INSTANCE,
		PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT | PIPE_REJECT_REMOTE_CLIENTS,
		1, 1024, 1024, NMPWAIT_USE_DEFAULT_WAIT, NULL
	);

	return pipeHandle;
}
/**
*  WaitForClient
*     Wait for the client to connect to a named pipe.
*    @param pipe - THe server pipe.
*    @return zero in the event of an error.
*/
int   WaitForClient(HANDLE pipe)
{
	return ConnectNamedPipe(pipe, NULL);
}
DWORD    WriteToPipe(HANDLE pipe, const void* pData, int nBytes)
{
	BOOL status;
	DWORD numBytesWritten;

	status = WriteFile(pipe, pData, nBytes, &numBytesWritten, NULL);
	if (status) {
		return numBytesWritten;
	}
	else {
		return 0;
	}



}
/**
* OpenClientPipe
*
*  Open a pipe as a client.  If the pipe does not yet exist, the function will retry
*  every 500ms.  Note that if you got the pipename wrong, this function will never terminate.
*  
*  @param basename - The base name of the pipe,  The full name will be 
*                     \\.\pipe\{basename} e.g. if basename is mypipe, the actual entity opened will
*                     be \\.\pipe\mypipe.
* @return HANDLE on failure this will be INVALID_HANDLE_VALUE  
*/
HANDLE OpenClientPipe(const char* basename)
{
	HANDLE h;
	wchar_t   fullName[MAX_PATH];
	swprintf(fullName, MAX_PATH, L"\\\\.\\pipe\\%s", basename);
	DWORD    err;
	do {
		h = CreateFile(
			fullName, GENERIC_READ, 0, NULL, OPEN_EXISTING, 0, NULL
			);
		err = GetLastError();
		if (h == INVALID_HANDLE_VALUE) {
			Sleep(500);
		}
	} while ((h == INVALID_HANDLE_VALUE) && ((err == ERROR_FILE_NOT_FOUND) || (err == ERROR_PATH_NOT_FOUND)));

	/* If control lands here, either we have a valid handle or there's an error that 
	   has nothing to do with the pipe not yet existing.
	  */

	return h;
}
/** 
 ReadFromPipe
    Read data froma pipe that has been opened.  Note that the functions in this module are 
	restricted to this being a client operation given the way the pipe is created by the server
	and how it is opened by the client.

	@param h      - Handle open on the pipe (from OpenClientPipe).
	@param pData  - Pointer to the buffer into which data will be read.
	@param nBytes - Size of the data buffer in bytes.
	@return       - Number of bytes read.
	@retval 0     - Error reading the pipe.
*/
DWORD    ReadFromPipe(HANDLE h, void* pData, int nBytes)
{
	DWORD bytesRead;

	if (ReadFile(h, pData, nBytes, &bytesRead, NULL)) {
		return bytesRead;
	} else {
		return 0;
	}

}
/**
*  ClosePipe
*    Called by client or server to close the handle on a named pipe.
*  @param h - pipe handle.. no longer valid for pipe operations after return.
*  @return int - nonzero for success.
*/
int ClosePipe(HANDLE h)
{
	return CloseHandle(h) ? 1 : 0;
}

#endif