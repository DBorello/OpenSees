/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.8 $
// $Date: 2006-11-30 23:51:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/Socket.h,v $
                                                                        
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains all the includes needed by files
// using the Berkeley socket stuff.

#ifndef Socket_h
#define Socket_h

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
  #include <winsock2.h>
#else
  #include <sys/socket.h>
  #include <sys/types.h>
  #include <netinet/in.h>
  #include <netinet/tcp.h>
  #include <arpa/inet.h>
  #include <netdb.h>
  #include <unistd.h>
  #include <strings.h>
#endif
}

#include <OPS_Stream.h>
#include <bool.h>

#define MAX_UDP_DATAGRAM 9126
#define MAX_INET_ADDR 28

#ifdef _WIN32
  typedef SOCKET socket_type;
  typedef int socklen_type;
  #define bzero(s,n) memset((s),0,(n))
  #define bcmp(s1,s2,n) memcmp((s1),(s2),(n))
#else
  typedef int socket_type;
  typedef socklen_t socklen_type;
#endif

extern int startup_sockets(void);
extern void cleanup_sockets(void);
extern int httpGet(const char *URL, const char *page, unsigned int port, char **dataPtr);

#endif 


