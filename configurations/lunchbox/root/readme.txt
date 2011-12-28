Copyright Dell Inc. 2004.  All rights reserved.	

Criticality: Recommended. Dell recommends applying this update during your next scheduled update cycle.  The update contains feature enhancements or changes that will help keep your system software current and compatible with other system modules (firmware, BIOS, drivers and software).

Title / Description:
Dell controller Linux CLI and SNMP agent

Contents:
1. Title / Description
2. Compatibility / Minimum Requirements
3. Release Highlights - Fixes
4. Installation
5. Known Limitations
6. History
7. General Information

Compatibility / Minimum Requirements:

WHQL certification:  	__ Yes	___ No		_X__ Not Applicable

Release Highlights - Features:
Initial Release

Release Highlights - Fixes:
Initial Release 

Installation:
-----------------------------------
Installing the CLI in Linux

1 - Log on as root or become superuser.
2 - If you are upgrading from an earlier version, check to see
if an afaapps kit already exists on your system, and if it does, 
remove it.
a. To determine if the afaapps kit is installed on your system, 
enter the following command:  	rpm -q afaapps
If the application software exists on your system, the rpm command
returns the following message, where VERSION_NUMBER is the actual
version of the application software:  afaapps-VERSION_NUMBER
b. If the application software already exists on your system, remove
it by entering the following command, replacing VERSION_NUMBER with
the version number returned by the rpm -q command:  
rpm -e afaapps-VERSION_NUMBER
3.  Download the afa-apps-snmp.2807420-a04.tar.gz package.
4.  Extract the package using "tar -xvzf <full path to the archive>"
This package will extract to the current directory and contains 
the following files:
afaapps-4.1-0.i386.rpm 
afasnmp-4.0-0.i386.rpm
readme.txt 
5.  As root, install the application package using "rpm -ivh 
afaapps-4.1-0.i386.rpm"


-----------------------------------
Installing the SNMP subagent in Linux

1.  Log on as root or become superuser.
2.  Check to see if the afasnmp SNMP Subagent for your controller
already exists on your system, and if it does, remove it.
a.  To determine if the afasnmp SNMP Subagent is installed on your
Linux system, enter the following command: rpm -q afasnmp
If the afasnmp SNMP Subagent exists on your Linux system, the "rpm"
command will return the following message, where "VERSION_NUMBER"
is the actual version of the SNMP software: afasnmp-"VERSION_NUMBER"
b.  If the afasnmp SNMP Subagent already exists on your Linux system,
remove it by entering the following command, replacing 
"VERSION_NUMBER" with the version number returned by the 
"rpm -q" command: rpm -e afasnmp-"VERSION_NUMBER"
3.  Install the afasnmp SNMP Subagent.
a.  Download the afa-apps-snmp.2807420-a04.tar.gz package.
b.  Extract the package using "tar -xvzf <full path to the archive>"
This package will extract to the current directory and contains 
the following files:
afaapps-4.1-0.i386.rpm
afasnmp-4.0-0.i386.rpm
readme.txt 
c.  As root, install the SNMP subagent using "rpm -ivh 
afasnmp-4.0-0.i386.rpm"

Known Limitations:
This utility is a legacy utility and there are some inherent limitations.
Disk set smart /logerr=1  (SCSI device)" command returns an error when used with the Crec Sata 6Ch controller.
Some enclosure commands fail on internal SCSI backplane and external enclosures.
Controller firmware save" is an old command not currently supported.
AFACLI generates an error message when delete a continer 32 or higher
When deleting a container, the last device from fdisk -l dissappears, and the system must be rebooted to have accurate status
"Enclosure set SCSIID" is not supported.











History: 
A02 Initial Release.                afa-apps-snmp.2.806076
A04 Minor bug fixes in afacli. afa-apps-snmp.2.807420



============================================================
General Information:

The Linux CLI allows a raid 1 to be created on a 2 disks configuration as well as a single disk configuration. Under a single disk configuration the array will not be redundant if a hard drive failure happens.

