

library(shiny)
library(sf)
library(leaflet)
library(ggplot2)
library(dplyr)
library(readr)
library(viridisLite)

# # -------------------------------------
# # Save as RData file (Run only once)
# # -------------------------------------
# drive_path <- "Z:/Projects/WP000023_Guinea_worm/Working/Som/data_NG/Shiny/"
# 
# adm1_shp_path <- file.path(drive_path, "geoBoundaries-NGA-ADM1-all", "geoBoundaries-NGA-ADM1.shp")
# adm1_csv_path <- file.path(drive_path, "NGA_adm1_estimates_dtp1_survey.csv")
# 
# # -----------------------------
# # Load shapefile (Admin 1 only)
# # -----------------------------
# cat("Loading Admin 1 shapefile...\n")
# shp_adm1 <- st_read(adm1_shp_path, quiet = TRUE)
# 
# # -----------------------------
# # Load CSV data (Admin 1 only)
# # -----------------------------
# cat("Loading Admin 1 CSV data...\n")
# adm1_data <- read_csv(adm1_csv_path, show_col_types = FALSE)


# -----------------------------
# Load pre-built RData (bundled with the app)
# -----------------------------
load("nga_dtp1_admin1.RData")
rdata_file <- "nga_dtp1_admin1.RData"

if (!file.exists(rdata_file)) {
  stop(paste("RData file not found:", rdata_file,
             "- make sure it is included in the deployment bundle."))
}

cat("Loading bundled RData...\n")
# load(rdata_file)   # brings shp_adm1 and adm1_data into the environment

# Sanity check
stopifnot(exists("shp_adm1"), exists("adm1_data"))

# -----------------------------
# Process ADM1 data -> %
# -----------------------------
process_adm_data <- function(df, upper_col, lower_col) {
  df <- df %>%
    mutate(
      Year = as.numeric(Year),
      ci_width = .data[[upper_col]] - .data[[lower_col]]
    ) %>%
    filter(Year >= 2000, Year <= 2024)
  
  if ("mean" %in% names(df) && all(df$mean >= 0 & df$mean <= 1, na.rm = TRUE)) {
    df$mean <- df$mean * 100
  }
  if ("sd" %in% names(df) && all(df$sd >= 0 & df$sd <= 1, na.rm = TRUE)) {
    df$sd <- df$sd * 100
  }
  if (all(df$ci_width >= 0 & df$ci_width <= 1, na.rm = TRUE)) {
    df$ci_width <- df$ci_width * 100
  }
  
  return(df)
}

adm1_data <- process_adm_data(adm1_data, "0.975quant", "0.025quant")

# Available years from the data
available_years <- sort(unique(adm1_data$Year))

# -----------------------------
# UI
# -----------------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }
      .well { background-color: #f8f9fa; border: 1px solid #dee2e6; }
      .form-group { margin-bottom: 15px; }
      .selectize-input { border-radius: 4px; border: 1px solid #ced4da; }
      .download-btn {
        background-color: #27ae60; color: white; border: none;
        padding: 10px 16px; border-radius: 4px; font-weight: bold;
        width: 100%; margin-top: 10px; cursor: pointer;
      }
      .download-btn:hover { background-color: #219150; color: white; }
    "))
  ),
  
  titlePanel(
    div(
      tags$i(class = "fas fa-syringe", style = "color: #e74c3c; margin-right: 10px;"),
      "Nigeria DTP1 Vaccination Coverage Estimates - Admin 1",
      style = "color: #2c3e50;"
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      div(class = "well",
          h4("Data Selection"),
          selectInput(
            "stat_type", "Estimate Type",
            choices = c("Mean (Coverage %)" = "mean",
                        "Standard Deviation (%)" = "sd",
                        "95% CI Width (%)" = "ci_width"),
            selected = "mean"
          ),
          selectInput(
            "year", "Year",
            choices = available_years,
            selected = max(available_years)
          ),
          downloadButton("download_map", "Download Map as PNG",
                         class = "download-btn")
      )
    ),
    
    mainPanel(
      width = 8,
      leafletOutput("map", height = 700)
    )
  )
)

# -----------------------------
# Server
# -----------------------------
server <- function(input, output, session) {
  
  stat_col <- reactive({
    switch(input$stat_type, "mean" = "mean", "sd" = "sd", "ci_width" = "ci_width")
  })
  
  stat_info <- reactive({
    switch(input$stat_type,
           "mean" = list(label = "Coverage", suffix = "%", palette = "viridis"),
           "sd" = list(label = "Standard Deviation", suffix = "%", palette = "cividis"),
           "ci_width" = list(label = "95% CI Width", suffix = "%", palette = "cividis"))
  })
  
  # Merged data for the current year
  merged_data <- reactive({
    req(input$year, input$stat_type)
    year_data <- adm1_data %>% filter(Year == input$year)
    if (nrow(year_data) == 0) {
      showNotification(paste("No data for year", input$year), type = "warning")
      return(NULL)
    }
    if (nrow(year_data) != nrow(shp_adm1)) {
      showNotification("Number of CSV rows does not match number of polygons", type = "error")
      return(NULL)
    }
    shp <- shp_adm1
    shp$value <- year_data[[stat_col()]]
    return(shp)
  })
  
  build_popup <- function(md, name_field, stat_label, suffix, year) {
    paste0(
      "<strong>Name:</strong> ", md[[name_field]], "<br>",
      "<strong>", stat_label, ":</strong> ", sprintf("%.3f", md$value), suffix, "<br>",
      "<strong>Year:</strong> ", year
    )
  }
  
  name_field <- reactive({
    md <- merged_data()
    if (is.null(md)) return("shapeName")
    if ("shapeName" %in% names(md)) "shapeName" else names(md)[1]
  })
  
  # Main interactive leaflet map
  output$map <- renderLeaflet({
    md <- merged_data()
    osm_url <- "https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png"
    osm_attribution <- '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
    
    if (is.null(md)) {
      return(leaflet() %>%
               addTiles(urlTemplate = osm_url, attribution = osm_attribution,
                        options = tileOptions(opacity = 0.4)) %>%
               setView(lng = 8, lat = 9, zoom = 6) %>%
               addControl("No data available", position = "topright"))
    }
    
    nf <- name_field()
    si <- stat_info()
    pal <- colorNumeric(si$palette, domain = md$value, na.color = "transparent")
    popup_content <- build_popup(md, nf, si$label, si$suffix, input$year)
    
    leaflet() %>%
      addTiles(urlTemplate = osm_url, attribution = osm_attribution,
               options = tileOptions(opacity = 0.4)) %>%
      addPolygons(
        data = md, fillColor = ~pal(value), fillOpacity = 0.85,
        color = "black", weight = 1, opacity = 0.8,
        popup = popup_content,
        highlightOptions = highlightOptions(weight = 3, color = "#666",
                                            fillOpacity = 0.9, bringToFront = TRUE),
        label = ~paste0(md[[nf]], ": ", sprintf("%.3f", md$value), si$suffix)
      ) %>%
      addLegend(position = "bottomright", pal = pal, values = md$value,
                title = paste0(si$label,
                               ifelse(si$suffix != "", paste0(" (", si$suffix, ")"), "")),
                opacity = 0.8, labFormat = labelFormat(suffix = si$suffix, digits = 3)) %>%
      setView(lng = 8, lat = 9, zoom = 6)
  })
  
  observe({
    md <- merged_data()
    if (!is.null(md)) {
      nf <- name_field()
      si <- stat_info()
      pal <- colorNumeric(si$palette, domain = md$value, na.color = "transparent")
      popup_content <- build_popup(md, nf, si$label, si$suffix, input$year)
      
      leafletProxy("map") %>%
        clearShapes() %>% clearControls() %>%
        addPolygons(
          data = md, fillColor = ~pal(value), fillOpacity = 0.85,
          color = "black", weight = 1, opacity = 0.8,
          popup = popup_content,
          highlightOptions = highlightOptions(weight = 3, color = "#666",
                                              fillOpacity = 0.9, bringToFront = TRUE),
          label = ~paste0(md[[nf]], ": ", sprintf("%.3f", md$value), si$suffix)
        ) %>%
        addLegend(position = "bottomright", pal = pal, values = md$value,
                  title = paste0(si$label,
                                 ifelse(si$suffix != "", paste0(" (", si$suffix, ")"), "")),
                  opacity = 0.8, labFormat = labelFormat(suffix = si$suffix, digits = 3))
    }
  })
  
  # ---- Download map as static PNG (no basemap, ggplot-based) ----
  output$download_map <- downloadHandler(
    filename = function() {
      paste0("DTP1_Admin1_", input$stat_type, "_", input$year, ".png")
    },
    content = function(file) {
      md <- merged_data()
      if (is.null(md)) {
        showNotification("No data available to download.", type = "error")
        return(NULL)
      }
      
      si <- stat_info()
      title_text  <- paste0("Nigeria DTP1 ", si$label, " (%) - ", input$year)
      legend_title <- paste0(si$label, " (%)")
      
      # Pick the matching viridis option letter
      # viridis = "D" (default), cividis = "C"
      viridis_option <- if (si$palette == "viridis") "D" else "C"
      
      # Build the static map (no basemap)
      p <- ggplot(md) +
        geom_sf(aes(fill = value), color = "black", linewidth = 0.3) +
        scale_fill_gradientn(
          colors = viridisLite::viridis(256, option = viridis_option),
          limits = c(min(md$value, na.rm = TRUE), max(md$value, na.rm = TRUE)),
          name = legend_title,
          labels = function(x) sprintf("%.1f%%", x),
          na.value = "transparent"
        ) +
        labs(title = title_text) +
        theme_void(base_size = 14) +
        theme(
          plot.title = element_text(
            hjust = 0.5, size = 18, face = "bold", color = "#2c3e50",
            margin = margin(b = 10)
          ),
          legend.position = "right",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          legend.key.height = unit(1.2, "cm"),
          plot.margin = margin(15, 15, 15, 15)
        )
      
      ggsave(
        filename = file,
        plot = p,
        width = 11,
        height = 8.5,
        dpi = 150,
        bg = "white"
      )
    }
  )
}

# -----------------------------
# Run App
# -----------------------------
shinyApp(ui, server)