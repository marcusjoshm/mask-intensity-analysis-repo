// ============================================================
// ROI_to_CellposeLabels.ijm
// 
// Creates a cellpose-style segmentation mask from ROI Manager.
// Each ROI is filled with a unique integer label (1, 2, 3, …).
// Background remains 0.
//
// Prerequisites:
//   1. A single image of the correct dimensions is open.
//   2. ROI Manager contains the cell-boundary ROIs.
//
// Output:
//   A 16-bit label image (supports up to 65,535 cells).
//   Save the result as .tiff to preserve integer labels.
// ============================================================

// --- Sanity checks -----------------------------------------------------------
if (nImages == 0)
    exit("No image is open. Open a .tiff of the correct dimensions first.");

nROIs = roiManager("count");
if (nROIs == 0)
    exit("ROI Manager is empty. Load your cell-boundary ROIs first.");

if (nROIs > 65535)
    exit("More than 65,535 ROIs — 16-bit label image cannot hold that many unique values.");

// --- Grab source image dimensions --------------------------------------------
title = getTitle();
w = getWidth();
h = getHeight();

// --- Create blank 16-bit label image -----------------------------------------
labelTitle = "CellposeLabels_" + title;
newImage(labelTitle, "16-bit black", w, h, 1);
selectWindow(labelTitle);

// --- Fill each ROI with its unique integer label -----------------------------
setBatchMode(true);          // speeds things up considerably

for (i = 0; i < nROIs; i++) {
    selectWindow(labelTitle);
    roiManager("select", i);

    // Label value = i + 1  (background stays 0)
    labelValue = i + 1;
    setColor(labelValue);
    fill();                  // fills inside the current ROI selection
}

// --- Clean up ----------------------------------------------------------------
run("Select None");
setBatchMode(false);

// --- Inform the user ---------------------------------------------------------
print("\\Clear");
print("=== Cell-Segmentation Label Mask Created ===");
print("Image       : " + labelTitle);
print("Dimensions  : " + w + " x " + h);
print("Total cells : " + nROIs);
print("Label range : 1 – " + nROIs);
print("");
print("Save as .tiff to preserve integer labels.");
print("(File > Save As > Tiff...)");
