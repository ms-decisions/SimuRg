## Author: Victor Sokolov
## First created: 2026-10-03
## Description: SimuRg
## Keywords: SimuRg tests
rm(list=ls())


#####--------------- Load functions and libraries ---------------#####
library(tidyverse)
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank(),
             legend.position = "top")


#####--------------- sg_translator tests ---------------#####
# rxode2 to Monolix
sg_translator("scripts/nlme/sg-translator/rxode_models/odes/model_PK_hv_iv.txt", to = "mlxtran",
              output_path = "test_model_PK_hv_iv_mlx.txt",
              dm_list = list(cmt = 1, adm = 1), output_vars = c("Cc_mgL", "Cc_mgL"))

sg_translator("scripts/nlme/sg-translator/rxode_models/macros/model_PK_2c.txt", to = "mlxtran",
              output_path = "test_model_PK_2c_mlx.txt", output_vars = c("Cc_mgL"))

sg_translator("scripts/nlme/sg-translator/rxode_models/macros/model_PK_2c.txt", to = "mlxtran",
              output_path = "test_model_PK_2c_mlx_nomacro.txt", output_vars = c("Cc_mgL"), macros = F)

sg_translator("scripts/nlme/sg-translator/rxode_models/macros/model_PK_2c.txt", to = "mlxtran",
              output_path = "test_model_PK_2c_mlx_nomacro_doses.txt", output_vars = c("Cc_mgL"), macros = F,
              dm_list = list(cmt = c(1, 1), adm = c(1, 2)))

sg_translator("scripts/nlme/sg-translator/rxode_models/macros/model_PK_2c.txt", to = "mlxtran",
              output_path = "test_model_PK_2c_mlx_doses.txt", output_vars = c("Cc_mgL"),
              dm_list = list(cmt = c(1, 1), adm = c(1, 2)))

# Monolix to rxode2
sg_translator("scripts/nlme/sg-translator/monolix_models/odes/model_PK_hv_iv.txt", to = "rxode",
              output_path = "test_model_PK_hv_iv_rx.txt")

### Validation
sg_translator("scripts/nlme/sg-translator/validation_monolix/danu_pk_1c.txt", to = "rxode",
              output_path = "val_danu_pk_1c_rx.txt")

sg_translator("scripts/nlme/sg-translator/validation_rxode/maritide_PK_model1.txt", to = "mlxtran",
              output_path = "val_maritide_PK_model1_mlx.txt", output_vars = c("Cc"))

sg_translator("scripts/nlme/sg-translator/validation_rxode/maritide_PK_model1.txt", to = "mlxtran",
              output_path = "val_maritide_PK_model1_mlx_nomacro.txt", output_vars = c("Cc"), macros = F)

sg_translator("scripts/nlme/sg-translator/validation_rxode/model-pcsk9-qsp.txt", to = "mlxtran",
              output_path = "val_model-pcsk9-qsp_mlx.txt", output_vars = c("PCSK9_CFB", "LDLc_CFB"))

sg_translator("scripts/nlme/sg-translator/validation_rxode/model-sglt2-qsp.txt", to = "mlxtran",
              output_path = "val_model-sglt2-qsp_mlx.txt", output_vars = c("UGE"),
              dm_list = list(cmt = c(1, 13), adm = c(1, 2)))


#####--------------- sg_converter tests ---------------#####
test_folder_2023 <- "scripts/nlme/2.1-sg-converter/monolix-2023/"
test_folder_2024 <- "scripts/nlme/2.1-sg-converter/monolix-2024/"
pro_name <- "pk-1c"

sg_converter(folder_path = test_folder_2023, proj_name = pro_name)
sg_converter(folder_path = test_folder_2024, proj_name = pro_name)
