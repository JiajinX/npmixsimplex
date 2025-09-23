library(readxl)

# Data clean for GDP compositions
gdpraw <- read_excel("data/raw/GDP.xlsx", skip = 2)
gdp2023 <- subset(gdpraw, select = c(Country, IndicatorName, `2023`))
gdp2023 <- clean_gdp(gdp2023)
# Income groups of countries
income <- read_excel("data/raw/OGHIST_2025_07_01.xlsx",
                     sheet = "Country Analytical History", skip = 5)
income <- income[6:223, c("Data for calendar year :", "2023")]
colnames(income) <- c("Area", "Income")
# Manually matching country names
income <- income %>%
  mutate(Area = recode(Area,
                       "Bahamas, The" = "Bahamas",
                       "Bolivia" = "Bolivia (Plurinational State of)",
                       "Hong Kong SAR, China" = "China, Hong Kong SAR",
                       "Macao SAR, China" = "China, Macao SAR",
                       "Yemen, Rep." = "Yemen",
                       "Venezuela, RB" = "Venezuela (Bolivarian Republic of)",
                       "Congo, Rep." = "Congo",
                       "Congo, Dem. Rep." = "D.R. of the Congo",
                       "Egypt, Arab Rep." = "Egypt",
                       "Korea, Dem. Rep." = "D.P.R. of Korea",
                       "Gambia, The" = "Gambia",
                       "Iran, Islamic Rep." = "Iran (Islamic Republic of)",
                       "Korea, Rep."  = "Republic of Korea",
                       "Kyrgyz Republic" = "Kyrgyzstan",
                       "Lao PDR" = "Lao People's DR",
                       "Micronesia, Fed. Sts." = "Micronesia (FS of)",
                       "Puerto Rico (U.S.)" = "Puerto Rico",
                       "Moldova" = "Republic of Moldova",
                       "St. Kitts and Nevis" = "Saint Kitts and Nevis",
                       "St. Lucia" = "Saint Lucia",
                       "Tanzania" = "U.R. of Tanzania: Mainland",
                       "São Tomé and Príncipe" = "Sao Tome and Principe",
                       "Slovak Republic" = "Slovakia",
                       "West Bank and Gaza" = "State of Palestine"))
gdp2023 <- dplyr::left_join(gdp2023, income)
gdp2023[gdp2023$Area == "Zanzibar", ]$Income <- "LM"
# saveRDS(gdp2023, "data/processed/gdp2023")

