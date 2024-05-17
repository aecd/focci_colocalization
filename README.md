# Focci colocalization

This is a FIJJ[1] macro script to automate focci colocalization in .HTD images from the ImageXpress spinning disk confocal.
The analysis pipeline follows the work described in Panichnantakul et al[2]. 

## How to use:
Install Comdet[3], and EzColocalization[4] according to the developer's instructions.

Simply, load the main.py into FIJI and run it. 
By including a background image in the same folder as the image, the script also performs a background subtraction (division).

## References:
1. Schindelin, J., Arganda-Carreras, I., Frise, E., Kaynig, V., Longair, M., Pietzsch, T., … Cardona, A. (2012) Fiji: an open-source platform for biological-image analysis. Nature Methods, 9(7), 676–682.
2. Panichnantakul et al (2021) DNA Repair 105, 103156
3. Comdet. https://github.com/UU-cellbiology/ComDet
4. Stauffer, Weston, Huanjie Sheng, and Han N. Lim. (2018) EzColocalization: An ImageJ plugin for visualizing and measuring colocalization in cells and organisms. Scientific reports 8.1: 15764.
