# Pulse Shape Discrimination (PSD) for the Musun experiment
This experiment measure the muon capture rate in deuterium to find a universal low energy constant for two nucleon weak interactions 
defined from EFT(effective field theory). Muon decay produce a large number of gamma rays and electrons The neutrons are detected by the 
ionization of the elastically scattered recoil protons following neutron, proton scattering whereas gamma rays are detected by the ionization 
of the recoil electrons following Compton scattering. Since heavier parti- cles have higher specific ionization, the protons produce more 
delayed fluorescence light compared to the electrons which results in a larger tail area for neutron pulses compared to gamma ray pulses This 
distinguishing feature is used for PSD.
The total area of the pulse is compared to the tail area to determine if the pulse was a neutron or a gamma ray. PSD ratio is the ratio of the 
tail area of the pulse to its total area. The figure below shows a plot of the PSD ratio versus the total area.</br>
<img width="429" alt="Screen Shot 2022-04-10 at 8 03 42 AM" src="https://user-images.githubusercontent.com/27436642/162617159-3d391da6-f693-4148-8beb-ce452e525246.png"></br>
This code makes template histograms by subtracting the pedestals and normalizing the sum of each pulse, to unit area.</br>
<img width="378" alt="Screen Shot 2022-04-10 at 8 09 37 AM" src="https://user-images.githubusercontent.com/27436642/162617359-bd589d76-e54f-4592-bf9f-7f74531f290d.png"></br>
Then each pulse is fitted using both templates with efficient and quick minimization called Brent Minimization technique. The minimization function 
is formed from the templates using the following formula:</br>
<img width="256" alt="Screen Shot 2022-04-10 at 8 13 00 AM" src="https://user-images.githubusercontent.com/27436642/162617475-b6427f91-f890-4275-9a90-743f820df263.png"></br>
The code uses MINUIT for fitting and based on the &chi;<sup>2</sup>/NDF we dertermine that a pulse is a neutron or a gamma pulse. 
Shown below is a neutron pulse fitted with a neutron template (left) and a gamma template (right):</br>
<img width="685" alt="Screen Shot 2022-04-10 at 8 16 44 AM" src="https://user-images.githubusercontent.com/27436642/162617636-26e03430-858d-42c9-a262-a676366c8931.png"></br>
This lower &chi;<sup>2</sup>/NDF with the neutron template fit clearly shows that it is a nuetron pulse.
## Instructions for compiling and running the code
#### Must have at C++ and ROOT 6.24/06  installed.
All code is in C and requires root to compile. Follow the simple instructions below</br>
& make clean</br>
& make</br>
& ./mu </br>
The basic code runs with the executable mu but this has a few options too (if you want less events).
