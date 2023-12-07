# dishSwSuite
Suite of software for the optical-geometric characterization of solar dishes

General information
	- Language: C++
	- Graphical User Interface: Qt
	- Library to make plot: Qwt
	- Library to manage images: OpenCV
	- Library to control GigE cameras: Aravis, Vimba

Warning
This software is offered as open source under the GNU General Public License v3.0
with the hope to be useful for someone. I apologize for the lack of accurate documentation, but at the present I have not time and funds for drafting that.
On the other hand, I am available to assist people interested on these software.

kFluxMapper
Software to control a GigE camera (by means of the Aravis library) for acquiring the images needed to the experimental evaluation of the profile of the solar radiation concentrated by the solar dish in its focal plane where a diffusive plane target is set.
Images of 
- Background
- oneSun (the radiation reflected by a plane mirror used to normalize the flux, obtaining the concentration factor)
- flux
can be shot and saved at one or more increasing exposure time

IWB – ImageWorkBench
This software processes the images acquired by kFluxMapper for composing the meta-image (double float) mapping the experimental flux.

kVISdish (master)
Software for infield experimental evaluation of 3D-shape of solar dishes according to the VIS approach: a suitably wide LCD monitor is placed orthogonally and centrally to the paraboloid axis at the point (optically) conjugated with the observation point, always placed along the paraboloid axis. The LCD monitor is connected to the “slave” PC, while the GigE camera, placed at the observation point, is connected to the neighbor “master” PC. KVISdish and kVISdishSlave have to be installed on the “master” and “slave” PC, respectively. “Master” and “slave” PCs are connected by Ethernet
KVISdish manages the workflow of the measurement:
1) send instruction to kVISdishSlave on the image that is to be shown on the LCD
2) controls the photo shooting by the GigE camera (Vimba library)
3) after the completion of the image acquisition, run the data processing

kVISdishSlave
Draws the image according to the information received by the “master”, and displays it on the LCD monitor.

SimulDish
Software for the evaluation of the flux distribution on the focal plane of a dish. The computing is based on the knowledge of the experimental 3D shape of the dish (obtained by the VISdish) as well as of the experimental Sun Conic Reflectance (SCR), i.e. the near-specular reflectance measured with a beam with divergence set to 4.7 mrad, which is the standard value for the Sun (half angle); SCR must be evaluate in the range of incidence and acceptance angle experienced in the dish.
The computing is dealt point by point in two steps: ray tracing of the central axis of the conic solar radiation hitting and reflected at the point P of the dish surface, until its crossing with the plane where is placed a virtual CCD; from this point the contribution to the flux of the point P is obtained by charging each virtual pixel with the amount of radiation described by SCR. Several different virtual CCD are considered at the same time in order to also determine the focal length. At the end, the flux contour-map and the value of the most relevant parameters (effective area, optimal focal length, intercept-factor, maximum and mean concentration,…) are displayed.
