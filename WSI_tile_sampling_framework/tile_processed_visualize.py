import os
import openslide
import h5py
from PIL import Image

import matplotlib.pyplot as plt
from PIL import Image, ImageDraw
import matplotlib as mpl
mpl.rcParams["figure.dpi"] = 300

import argparse

parser = argparse.ArgumentParser()
#parser.add_argument('--organ', type= str , default = 'ACC', help = 'WSI data directory folder')

args = parser.parse_args()

def main() : 

    #root = '/mnt/fileserver1_data/nfs/shared/Pathology/Prostate_Busan/busan_set8/'
    #names = ['DB-000508-B3']
    #names = [word[:-3] for word in os.listdir('/mnt/fileserver1_data/nfs/shared/Pathology/Prostate_Busan/embedded_features/coords')]
    
    root = '/mnt/fileserver1_data/nfs/shared/Pathology/Prostate_Busan/'
    

    names = []
    states = ["busan_set{}".format(i) for i in range(1, 16)]
    for state in states :
        for slides in os.listdir(os.path.join(root, state)) : 
            if slides.endswith(".svs") : 
                try : 
                    if slides.split('-')[2][0] == 'A' : ### 일부 sampling (A의 말초대와 중심부인 B~D만 sampling, E,F와 기타 부위들 전부 샘플링 x)
                        final_directory = os.path.join(root,state,slides)
                        names.append(final_directory)
                except : 
                    continue

    
    #names = ['/mnt/fileserver1_data/nfs/shared/Pathology/Prostate_Busan/busan_set8/DB-000524-F6.svs']
    #names = ['/mnt/fileserver1_data/nfs/shared/Pathology/Prostate_Busan/busan_set1/DB-000020-A1.svs']
    #names = ['/mnt/fileserver1_data/nfs/shared/Pathology/Prostate_Busan/busan_set10/DB-000660-B4.svs']
    
    thumbnail_level = -1

    for name in names : 
        print("VISUALUIZATION {} ".format(name))
        slide_name = name.split('/')[-1]
        slide_path = name
        name = slide_name[:-4]
        
        #slide_name = name + '.svs'
        h5_path = '/mnt/fileserver1_data/nfs/shared/Pathology/Prostate_Busan/embedded_features/coords_all/{}.h5'.format(name)
        #h5_path = '/home/yscho/tissue_cell_segmentation_PNU/Hover_next/DB-000524-F6/coords/DB-000524-F6.h5'
        #h5_path = '/home/yscho/tissue_cell_segmentation_PNU/Hover_next/DB-000660-B4/coords/DB-000660-B4.h5'
        
        with h5py.File(h5_path, "r") as h5_file:
            tile_coords = h5_file['tile_coords'][:]
        
        print("Slide Path", os.path.join(root, slide_name))
        #slide = openslide.OpenSlide(os.path.join(root, slide_name))
        slide = openslide.OpenSlide(slide_path)
        
        print("Original Size", slide.level_dimensions)
        if len(slide.level_dimensions) > 1:
            thumbnail_size = slide.level_dimensions[thumbnail_level] 
        else : 
            thumbnail_size = (slide.level_dimensions[0][0]//32,slide.level_dimensions[0][1]//32)
        print("Thumbnail Size", thumbnail_size)
        thumbnail = slide.get_thumbnail(thumbnail_size).convert("RGB")  # RGB로 변환

        # 2. 원본 이미지 크기
        w_orig, h_orig = slide.dimensions
        print(f"Original size: ({w_orig}, {h_orig})")
        print(f"Thumbnail size: {thumbnail_size}")

        # 3. 축소 비율 계산
        scale_x = thumbnail_size[0] / w_orig
        scale_y = thumbnail_size[1] / h_orig
        print(f"Scale factors: x={scale_x}, y={scale_y}")

        # 4. ImageDraw 객체를 루프 밖에서 한 번만 생성
        draw = ImageDraw.Draw(thumbnail)

        # 5. 여러 좌표에 대해 BBox 그리기
        for coord in tile_coords: 
            # 원본 좌표
            x_orig, y_orig = coord

            # 썸네일 좌표로 변환
            x_thumb = int(x_orig * scale_x)
            y_thumb = int(y_orig * scale_y)

            # BBox 크기 조정
            bbox_size = 512  # 원본 기준 크기
            bbox_width = int(bbox_size * scale_x)
            bbox_height = int(bbox_size * scale_y)

            # 좌표 및 BBox 범위 체크
            if (x_thumb < 0 or y_thumb < 0 or 
                x_thumb + bbox_width > thumbnail_size[0] or 
                y_thumb + bbox_height > thumbnail_size[1]):
                pass
            else:
                # 빨간색 사각형 그리기
                draw.rectangle(
                    [(x_thumb, y_thumb), (x_thumb + bbox_width, y_thumb + bbox_height)],
                    outline="red", width=1
                )

        # 6. 시각화 (루프 밖에서 한 번만 호출)
        plt.imshow(thumbnail)
        plt.axis('off')
        #plt.title("Thumbnail with Multiple BBoxes")
        #plt.show()
        plt.savefig('/mnt/fileserver1_data/nfs/shared/Pathology/Prostate_Busan/embedded_features/thumbnail_patch_all/{}.png'.format(name))
        #plt.savefig('/home/yscho/tissue_cell_segmentation_PNU/Hover_next/DB-000660-B4/thumbnail_patch/DB-000660-B4.png')
        
if __name__ == "__main__" :
    main()
