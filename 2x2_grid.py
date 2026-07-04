import os
from PIL import Image, ImageChops
import matplotlib.pyplot as plt


my_plots = [
    "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data_for_publication/plots/mybaits_sample4_coverage_map.png",
    "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data_for_publication/plots/mybaits_sample18_coverage_map.png",
    "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data_for_publication/plots/swift_sample4_coverage_map.png",
    "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data_for_publication/plots/swift_sample18_coverage_map.png",
]

final_output = "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data_for_publication/plots/sample4_18_2x2__final_combined_publication_figure.png"


def trim_white_borders(image_path):
    """Automatically crops out excess white margins around a plot."""
    img = Image.open(image_path)
    # Convert to RGB if it's in RGBA mode
    if img.mode == 'RGBA':
        background = Image.new('RGB', img.size, (255, 255, 255))
        background.paste(img, mask=img.split()[3])
        img = background

    # Find the bounding box of non-white areas
    bg = Image.new(img.mode, img.size, (255, 255, 255))
    diff = ImageChops.difference(img, bg)
    bbox = diff.getbbox()

    return img.crop(bbox) if bbox else img


def create_2x2_panel(image_paths, output_path, panel_labels=['A', 'B', 'C', 'D']):
    """Combines 4 images into a 2x2 grid panel with clean spacing."""
    if len(image_paths) != 4:
        raise ValueError("Exactly 4 image paths must be provided.")

    # Trim and load images
    images = [trim_white_borders(path) for path in image_paths]

    # Create the figure grid
    fig, axes = plt.subplots(2, 2, figsize=(16, 8.5), layout="constrained")
    axes_flat = axes.flatten()

    for i, ax in enumerate(axes_flat):
        # Display image in its respective grid cell
        ax.imshow(images[i])
        ax.axis('off')  # Hide axis ticks/frames

        # Add publication-ready panel labels (A, B, C, D) in the top-left corner
        ax.text(
            -0.02, 1.05, panel_labels[i],
            transform=ax.transAxes,
            fontsize=20,
            fontweight='bold',
            va='top',
            ha='right'
        )

    # Minimize spacing between the subplots
    plt.tight_layout(pad=1.0, w_pad=2.0, h_pad=2.0)

    # Save the combined panel with a high DPI
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Publication panel successfully saved to: {output_path}")


# ==========================================
# HOW TO RUN IT:
# ==========================================
# Define the paths to your 4 previously generated PNGs

# Run the merger
if all(os.path.exists(p) for p in my_plots):
    create_2x2_panel(my_plots, final_output)
else:
    print("Error: One or more of the 4 source image paths do not exist yet.")