这是对Tomographic Light Field Synthesis-Code的注释

所有的代码来自于：

* Layered 3D: Tomographic Image Synthesis for Attenuation-based Light Field and High Dynamic Range Displays
* Source Code and Light Field Datasets 
* Gordon Wetzstein
* Siggraph 2011

我发现要想比较清楚地弄明白光场的计算过程，有必要把这些代码完全读懂，而读懂代码最好的方法就是全文注释一遍。

阅读顺序：

1. reconstructLayers.m
2. precomputeSparsePropagationMatrixLayers3D.m

硬件的设置参考网页：

http://displayblocks.org/diycompressivedisplays/tensordisplays/

论文的原文在：

http://alumni.media.mit.edu/~dlanman/research/Layered3D/

光场显示的简介

用常规的显示器，例如液晶或者投影仪，每个像素所发出的光是弥散光，指向各个方向。在观察这样的显示器时，只能看到二维平面的图像。即使是使用了左右眼具有视差的技术，例如主动快门或者红蓝眼镜之类，使左右看到不同的图像，也只是给人以一定的“立体感”，而非真正产生了立体像。

光场显示技术用是要使显示器像素不仅仅表现强度的变化，还要使每个像素带有光线方向的信息。于是人眼在观察的时候，还可以获得光学上的立体效果。例如可以使用裸眼观察，或者当眼睛调节焦距时能够有远近不同时物体呈现模糊或清晰的区别。