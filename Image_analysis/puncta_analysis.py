import os
import cv2
import numpy as np
import pandas as pd
from skimage.measure import label as skimage_label, regionprops
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from skimage.filters import threshold_otsu

def otsu_thresholding(image, scale=1.0):
    """Apply scaled Otsu thresholding to an image to identify puncta."""
    threshold = threshold_otsu(image)
    scaled_threshold = threshold * scale
    binary_image = image > scaled_threshold
    return binary_image

def show_image(title, image):
    plt.figure(figsize=(10, 10))
    plt.title(title)
    plt.imshow(image, cmap='gray')
    plt.show()

def find_puncta(image, replicate,condition,image_folder,cell_name,  min_size=40.15):
    """Identify puncta in the thresholded image and return their average size."""
    labeled_array, num_features = skimage_label(image, return_num=True)
    regions = regionprops(labeled_array)
    filtered_areas = [region.area for region in regions if region.area > min_size]

    filtered_array = np.isin(labeled_array, [region.label for region in regions if region.area > min_size]) * labeled_array
    filtered_num_features = len(filtered_areas)
    average_size = np.mean(filtered_areas) if filtered_areas else 0

    # Physical conversion factors
    pixel_size_um = 78.38 / 2221  # microns per pixel
    pixel_area_um2 = pixel_size_um ** 2  # square microns per pixel

    # Convert pixel areas to square microns
    filtered_areas_um2 = [area * pixel_area_um2 for area in filtered_areas]
    average_size_um2 = np.mean(filtered_areas_um2) if filtered_areas_um2 else 0

    return filtered_array, filtered_num_features, average_size_um2, filtered_areas_um2

def calculate_overlap(atg9_labeled, ha_labeled,overlap_sizes_df, replicate,condition,image_folder,cell_name):
    """Calculate the number and adjusted average size of HA puncta overlapping with ATG9 puncta."""

    # ha_labeled is an array of 0,1 --> connected 1s are puncta
    ha_regions = regionprops(ha_labeled)
    atg9_regions = regionprops(atg9_labeled)
    # ha_labeled is the binary version of our image, already tresholded
    # regionprops picks out the regions of interest --> list of regions that are connected (region1, region2,...)

    overlap_count = 0
    overlap_areas = []

    pixel_size_um = 78.38 / 2221  # microns per pixel
    pixel_area_um2 = pixel_size_um ** 2  # square microns per pixel

    for ha_region in ha_regions:

        # iterates through the selected areas; each region has a number, here we are selecting for it
        ha_mask = ha_labeled == ha_region.label #ha_region.label --> is a 1,2,3,4, points to the necessary region we want to select
        overlap_mask = np.logical_and(ha_mask, atg9_labeled > 0)

        # select only the region of interest (select only are of ha_region)

        if np.any(overlap_mask):
            overlap_count += 1

            total_overlap_area = 0
            # iterate through all atg9 regions to find the ones that overlap, if they overlap this is appened on ATG9_overlap_area
            ha_min_row, ha_min_col, ha_max_row, ha_max_col = ha_region.bbox

            # Iterate through ATG9 regions that intersect with the HA region bounding box
            for atg9_region in atg9_regions:
                atg9_min_row, atg9_min_col, atg9_max_row, atg9_max_col = atg9_region.bbox
                if (atg9_min_row < ha_max_row and atg9_max_row > ha_min_row and
                    atg9_min_col < ha_max_col and atg9_max_col > ha_min_col):

                    atg9_mask = atg9_labeled == atg9_region.label
                    intersection_mask = np.logical_and(ha_mask, atg9_mask)
                    overlap_area = np.sum(intersection_mask)

                    if overlap_area > 0:
                        # adding all the areas that intersect with this HA_puncta (summing together)
                        total_overlap_area += overlap_area
            # one large array with all the ares --> will be averages

            if total_overlap_area>40.15: # min size
                overlap_areas.append(total_overlap_area)

                overlap_sizes_df.extend([{'Replicate': replicate, 'Condition': condition, 'Image': image_folder,
                                      'Cell': cell_name, 'Puncta_Size': total_overlap_area * pixel_area_um2}])

    average_overlap_size = np.mean(overlap_areas) if overlap_areas else 0
    average_overlap_size = average_overlap_size * pixel_area_um2

    return overlap_count, average_overlap_size, overlap_sizes_df

def calculate_pearson(image1, image2):
    """Calculate the Pearson correlation coefficient between two images."""
    flat_image1 = image1.flatten()
    flat_image2 = image2.flatten()

    if np.all(flat_image1 == flat_image1[0]) or np.all(flat_image2 == flat_image2[0]):
        null_pearson = 0
        return null_pearson

    pearson_coef, _ = pearsonr(flat_image1, flat_image2)
    return pearson_coef

def process_images(folder_path, threshold_scale_atg9=1.0, threshold_scale_ha=1.0, show_images=True, min_particle_size=40.15):
    threshold_scale_atg9 = max(0.1, min(threshold_scale_atg9, 5.0))
    threshold_scale_ha = max(0.1, min(threshold_scale_ha, 5.0))

    data = []
    puncta_sizes_ATG9 = []
    puncta_sizes_HA = []

    overlap_sizes_df = []

    for replicate in os.listdir(folder_path):
        print(replicate)
        replicate_path = os.path.join(folder_path, replicate)
        if os.path.isdir(replicate_path):
            for condition in os.listdir(replicate_path):
                print(condition)
                condition_path = os.path.join(replicate_path, condition)
                if os.path.isdir(condition_path):
                    for image_folder in os.listdir(condition_path):
                        image_path = os.path.join(condition_path, image_folder)
                        if os.path.isdir(image_path):
                            for cell_folder in os.listdir(image_path):
                                cell_path = os.path.join(image_path, cell_folder)
                                if os.path.isdir(cell_path):
                                    atg9_image = ha_image = None
                                    atg9_thresh = ha_thresh = None
                                    atg9_count = ha_count = 0
                                    atg9_avg_size = ha_avg_size = 0
                                    cell_name = cell_folder

                                    for file_name in os.listdir(cell_path):
                                        if file_name.endswith(".tif"):
                                            if "ATG9" in file_name:
                                                channel = "ATG9"
                                            elif "HA" in file_name:
                                                channel = "HA"
                                            else:
                                                print(f"Skipping file {file_name} with unexpected naming format.")
                                                continue

                                            image_path_full = os.path.join(cell_path, file_name)
                                            print(f"Processing file: {file_name} for channel: {channel}")

                                            if channel == "ATG9":
                                                atg9_image = cv2.imread(image_path_full, cv2.IMREAD_GRAYSCALE)
                                                if atg9_image is not None:
                                                    atg9_thresh = otsu_thresholding(atg9_image, scale=threshold_scale_atg9)
                                                    atg9_puncta, atg9_count, atg9_avg_size, filtered_areas_ATG9 = find_puncta(atg9_thresh, replicate,
                                                                                                                           condition, image_folder, cell_name, min_size=40.15)

                                                    puncta_sizes_ATG9.extend([{'Replicate': replicate,
                                                                                'Condition': condition,
                                                                                'Image': image_folder,
                                                                                'Cell': cell_name,
                                                                                'Puncta_Size': area} for
                                                                               area in filtered_areas_ATG9])

                                                    if show_images:
                                                        show_image(f"ATG9 Threshold - {file_name}", atg9_thresh)
                                                    print(f"ATG9 puncta found: {atg9_count}")
                                                else:
                                                    print(f"Failed to read {image_path_full}")
                                            elif channel == "HA":
                                                ha_image = cv2.imread(image_path_full, cv2.IMREAD_GRAYSCALE)
                                                if ha_image is not None:
                                                    ha_thresh = otsu_thresholding(ha_image, scale=threshold_scale_ha)
                                                    ha_puncta, ha_count, ha_avg_size, filtered_areas_HA = find_puncta(ha_thresh, replicate,
                                                                                                                    condition, image_folder, cell_name, min_size=40.15)

                                                    puncta_sizes_HA.extend([{'Replicate': replicate,
                                                                               'Condition': condition,
                                                                               'Image': image_folder,
                                                                               'Cell': cell_name,
                                                                               'Puncta_Size': area} for
                                                                              area in filtered_areas_HA])

                                                    if show_images:
                                                        show_image(f"HA Threshold - {file_name}", ha_thresh)
                                                    print(f"HA puncta found: {ha_count}")
                                                else:
                                                    print(f"Failed to read {image_path_full}")

                                    if atg9_image is not None and ha_image is not None:
                                        pearson_coef = calculate_pearson(atg9_thresh, ha_thresh)
                                    else:
                                        pearson_coef = None
                                        print(f"Pearson calculation skipped due to missing data for cell {cell_name} in image folder {image_folder}")

                                    if atg9_thresh is not None and ha_thresh is not None:
                                        overlap_count, avg_overlap_size,overlap_sizes_df = calculate_overlap(atg9_puncta, ha_puncta, overlap_sizes_df,
                                                                                                             replicate, condition, image_folder, cell_name)
                                    else:
                                        print(f"Missing data for cell {cell_name} in image folder {image_folder}")
                                        overlap_count = 0
                                        avg_overlap_size = 0

                                    data.append({
                                        'Replicate': replicate,
                                        'Condition': condition,
                                        'Image': image_folder,
                                        'Cell': cell_name,
                                        'ATG9_Puncta': atg9_count,
                                        'HA_Puncta': ha_count,
                                        'Overlap_Puncta': overlap_count,
                                        'ATG9_Avg_Size': atg9_avg_size,
                                        'HA_Avg_Size': ha_avg_size,
                                        'Overlap_Avg_Size': avg_overlap_size,
                                        'Pearson_Coefficient': pearson_coef
                                    })

    df = pd.DataFrame(data)
    puncta_sizes_ATG9 = pd.DataFrame(puncta_sizes_ATG9)
    puncta_sizes_HA=pd.DataFrame(puncta_sizes_HA)

    overlap_sizes_df = pd.DataFrame(overlap_sizes_df)
    return df, puncta_sizes_ATG9,puncta_sizes_HA, overlap_sizes_df

def main():
    folder_path = ''
    threshold_scale_atg9 = 2.5  # Adjust this value to change the threshold scale for ATG9 channel
    threshold_scale_ha = 3.8  # Adjust this value to change the threshold scale for HA channel
    show_images = False  # Set to True to display images
    min_particle_size = 40.15  # Set the minimum particle size threshold

    df, puncta_sizes_ATG9,puncta_sizes_HA, overlap_sizes_df = process_images(folder_path, threshold_scale_atg9, threshold_scale_ha, show_images, min_particle_size)
    print(df)
    df.to_csv('puncta_quantification_results_ATG9_HA.csv', index=False)
    print("Data saved to puncta_quantification_results_ATG9_HA.csv")

    puncta_sizes_ATG9.to_csv('all_puncta_sizes_ATG9.csv', index=False)
    puncta_sizes_HA.to_csv('all_puncta_sizes_HA.csv', index=False)

    print("Data saved to all_puncta_sizes_ATG9_HA.csv")

    overlap_sizes_df.to_csv("overlap_sizes_df_ATG9_HA.csv", index=False)

if __name__ == "__main__":
    main()
