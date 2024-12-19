#!/bin/bash

set -e
for dir in Logs PDFs Images; do
  if [ -d "$dir" ]; then
    echo "Directory $dir already exists, skipping creation."
  else
     mkdir "$dir"
     chmod -R ug=rwx,o=rx "$dir"
     echo "Directory $dir created successfully."
  fi
done

run_workbook_template() {
  local template_name=$1
  local template_location=$2
  local template_count=$3
  echo "############ ${template_location}/${template_count} ${template_name}.R Starting... ########"
  Rscript ${template_name}.R 2>&1 | tee ./Logs/${template_name}.log; exit_status=${PIPESTATUS[0]}
  if [ $exit_status -ne 0 ]; then
    exit $exit_status
  fi
  if [ -f Rplots.pdf ]; then
    mv Rplots.pdf ./PDFs/${template_name}.pdf
  fi
  if [ -f Plot_PCA_3D.html ]; then
    mv Plot_PCA_3D.html ./Images/${template_name}_Plot_PCA_3D.html
  fi
  find . -maxdepth 1 -type f \( -iname '*.png' -o -iname '*.jpg' \) -exec mv {} Images/ \;
  echo "############ ${template_location}/${template_count} ${template_name}.R Completed! ########"
}
############ 1/24 template_ccbr1223_human_CleanRawCounts.R ########
run_workbook_template "template_ccbr1223_human_CleanRawCounts" "1" "24"

############ 2/24 template_ccbr1223_human_FilteredCounts.R ########
run_workbook_template "template_ccbr1223_human_FilteredCounts" "2" "24"

############ 3/24 template_ccbr1223_human_NormalizedCounts.R ########
run_workbook_template "template_ccbr1223_human_NormalizedCounts" "3" "24"

############ 4/24 template_ccbr1223_mouse_CleanRawCounts.R ########
run_workbook_template "template_ccbr1223_mouse_CleanRawCounts" "4" "24"

############ 5/24 template_ccbr1223_mouse_FilteredCounts.R ########
run_workbook_template "template_ccbr1223_mouse_FilteredCounts" "5" "24"

############ 6/24 template_ccbr1223_mouse_NormalizedCounts.R ########
run_workbook_template "template_ccbr1223_mouse_NormalizedCounts" "6" "24"

############ 7/24 template_ccbr1223_human_BatchCorrectedCounts.R ########
run_workbook_template "template_ccbr1223_human_BatchCorrectedCounts" "7" "24"

############ 8/24 template_ccbr1223_human_DEG.R ########
run_workbook_template "template_ccbr1223_human_DEG" "8" "24"

############ 9/24 template_ccbr1223_human_DEGGeneList.R ########
run_workbook_template "template_ccbr1223_human_DEGGeneList" "9" "24"

############ 10/24 template_ccbr1223_human_ExpressionHeatmap.R ########
run_workbook_template "template_ccbr1223_human_ExpressionHeatmap" "10" "24"

############ 11/24 template_ccbr1223_human_GSEAPreranked.R ########
run_workbook_template "template_ccbr1223_human_GSEAPreranked" "11" "24"

############ 12/24 template_ccbr1223_human_PCA3D.R ########
run_workbook_template "template_ccbr1223_human_PCA3D" "12" "24"

############ 13/24 template_ccbr1223_human_Volcano.R ########
run_workbook_template "template_ccbr1223_human_Volcano" "13" "24"

############ 14/24 template_ccbr1223_mouse_BatchCorrectedCounts.R ########
run_workbook_template "template_ccbr1223_mouse_BatchCorrectedCounts" "14" "24"

############ 15/24 template_ccbr1223_mouse_DEG.R ########
run_workbook_template "template_ccbr1223_mouse_DEG" "15" "24"

############ 16/24 template_ccbr1223_mouse_DEGGeneList.R ########
run_workbook_template "template_ccbr1223_mouse_DEGGeneList" "16" "24"

############ 17/24 template_ccbr1223_mouse_ExpressionHeatmap.R ########
run_workbook_template "template_ccbr1223_mouse_ExpressionHeatmap" "17" "24"

############ 18/24 template_ccbr1223_mouse_GSEAPreranked.R ########
run_workbook_template "template_ccbr1223_mouse_GSEAPreranked" "18" "24"

############ 19/24 template_ccbr1223_mouse_PCA3D.R ########
run_workbook_template "template_ccbr1223_mouse_PCA3D" "19" "24"

############ 20/24 template_ccbr1223_mouse_Volcano.R ########
run_workbook_template "template_ccbr1223_mouse_Volcano" "20" "24"

############ 21/24 template_ccbr1223_human_GSEAFiltered.R ########
run_workbook_template "template_ccbr1223_human_GSEAFiltered" "21" "24"

############ 22/24 template_ccbr1223_human_GSEAVisualization.R ########
run_workbook_template "template_ccbr1223_human_GSEAVisualization" "22" "24"

############ 23/24 template_ccbr1223_mouse_GSEAFiltered.R ########
run_workbook_template "template_ccbr1223_mouse_GSEAFiltered" "23" "24"

############ 24/24 template_ccbr1223_mouse_GSEAVisualization.R ########
run_workbook_template "template_ccbr1223_mouse_GSEAVisualization" "24" "24"

