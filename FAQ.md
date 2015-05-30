---
layout: page
title: FAQ
---

### General:

**Q** Which computer labs on campus have Virtual Box?
>**A:** The following rooms have Virtual Box installed:  

>* 2060 SciLab
* 2101 SCC
* 182 Shields  

>The SciLab rooms have USB 3, but they also have classes. You can see when all those rooms are free [here](http://computerrooms.ucdavis.edu/rooms/available/) (click on "Next Day" to see the next 7 days).


**Q:** When is my next assignment due?
>**A:** Assignments are due by **9AM** on the date indicated on the [Assignments page](http://jnmaloof.github.io/BIS180L_web/assignments/).


**Q:** I can't run my virtual machine on my personal computer. Help!  
>**A:** Make sure you have Virtual Box downloaded and installed. Downloads for various files systems can be found [here](https://www.virtualbox.org/wiki/Downloads). Then try opening the `.vbox` file again.



**Q:** Virtual Box has *fill in frustrating error message here*. What can I do?
>**A:** Try the following: listed in order of severity.

>1. Check system architecture.
	+ Open Ubuntu Settings while the machine is not running
	+ In the General tab, it should say Version Ubuntu (64-bit).
	+ If it doesn't, change it.
2. Check memory settings
	+ Your machine could be running out of memory
	+ In Settings > System tab, change base memory to say 2500 MB.
3. If neither of these work
	+ Reboot the Virtual Box, holding down shift as soon as reboot starts.  
	+ Choose Advanced Options for Ubuntu. 
	+ Choose Ubuntu, with Linux 3.1.16.0-33-generic (recovery mode).  
	+ Select fsck. 
	+ Reboot again.
