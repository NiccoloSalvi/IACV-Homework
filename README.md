# 3D Geometry Recovery of a Rectangular Parallelepiped

This project was developed for the *Image Analysis and Computer Vision* (IACV) course – A.Y. 2024/2025.
It addresses the recovery of 3D geometric properties from a single image of a parallelepiped-shaped piece of furniture, using both theoretical and computational tools in MATLAB.

---

## 📄 Assignment Overview

Given:

* A **single image** of a parallelepiped (furniture) taken with an **uncalibrated, zero-skew camera**
* The **width** of the parallelepiped (`l = 1`) along the X-axis
* Unknown **depth** (`m`) and **height** (`h`)
* A visible **horizontal circumference** (on the X-Y plane)
* A visible **horizontal curve** located at `z = h/2`

Tasks:

1. Extract straight lines, the image of the circle `C`, and the curve `S`
2. Estimate the **vanishing line** and perform **rectification** of the horizontal plane
3. Calibrate the camera and compute the **depth** and **height** of the object
4. Rectify and extract **X-Y coordinates** from the curve `S`
5. Localize the **camera** with respect to the object
6. Visualize the **3D reconstructed model**

---

## 🛠️ Implementation

The project was implemented entirely in **MATLAB**, using:

* Feature extraction: `Canny` edge detector and `Hough` transform
* Geometric reasoning and projective geometry for vanishing points and rectification
* Direct Linear Transformation (DLT) for camera localization
* Visualization of results through rectified images and 3D rendering

---

## 📁 Repository Structure

```
📆 root/
├── F1/               # Extraction of straight lines
├── G1–G6/            # Geometry tasks (vanishing points, calibration, 3D model, etc.)
├── data/             # Input data, line coordinates, calibration data
├── images/           # Input and output images (rectified, labeled, etc.)
├── model3D/          # MATLAB scripts to build and render the 3D model
├── assignment.pdf    # Official project description
├── report.pdf        # Detailed report of the implementation and results
└── LICENSE           # License file
```

## 📌 References

* 📄 [`assignment.pdf`](assignment.pdf) – Project requirements
* 📝 [`report.pdf`](report.pdf) – Full implementation and results
