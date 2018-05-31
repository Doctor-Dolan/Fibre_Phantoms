function [] = Dice;

MASK=load_nii('Density.nii');
MAP=load_nii('Ushaped_0_VolumeFractions.nii');

mask=MASK.img;          %Loads as 16int
mask=sum(mask,4);       %Converts to double

map=MAP.img;
map=map(:,:,:,2);         %Select correct fibre from VF slices

mask(mask>0)=1;
map(map>0)=1;           %Data = 1 where data != 0

tmp=mask.*map;          %Dice Numerator
sumINT=sum(tmp(:));

sumMask=sum(mask(:));   %Dice Denominators
sumMap=sum(map(:));

Dice=2*sumINT/(sumMap+sumMask)
