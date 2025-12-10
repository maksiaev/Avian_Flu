# install.packages("tidyr")

# Load required libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
data <- read_excel("C:/Users/maksi/Documents/Statistics/Projects/Avian_Flu_Files/Cats/Cat_Table2_R.xlsx")


# Preview the data
print("Data summary:")
print(head(data))
print(paste("Total rows:", nrow(data)))

# Remove rows with NA in week ending date
data <- data %>% 
  filter(complete_collection_date == 1, exposure != "avian")

# For stacked lollipops, we need to count occurrences per week/genotype/exposure combination
# and create a position variable for stacking
data_stacked <- data %>%
  arrange(collection_date, exposure, genotype) %>%
  group_by(collection_date, exposure) %>%
  mutate(position = row_number()) %>%
  ungroup()

print(data_stacked)

#cols<- c("B3.13"= "orange", "B3.13 (inferred)"= "#FED8B1", "D1.1"= "#7F00FF", "D1.1 (inferred)"= "#D6B4FC",
 #        "B3.2"= "#895129", "D1.3"= "steelblue", "unknown"= "grey")

cols2<- c("B3.13"= "#2069f0", "D1.1"= "#9100bd", "B3.7" = "#219c8c", "B3.6"="#4cfe78", "B3.2"="#fed976",
          "A3"= "#ff7b00", "B3.5"="#f57a93")

start <- as.Date("2021-11-01", "%Y-%m-%d")
end <- as.Date("2025-08-01", "%Y-%m-%d")

# data_stacked$genotype <- factor(data_stacked$genotype, levels = c(""))

excel_origin <- "1899-12-30"

# Function to find the next desired end-of-week date (e.g., the first Saturday after a start date)
get_eow_dates <- function(start_date, end_date) { #}, day_of_week = "Saturday") {
  # Convert start and end dates to Date objects if they aren't already
  start_date <- as.Date(start_date, "%Y-%m-%d")
  end_date <- as.Date(end_date, "%Y-%m-%d")
  
  # Find the first end-of-week day (e.g., Saturday) on or after the start date
  first_eow <- as.Date(start_date + (which(weekdays(seq.Date(start_date, start_date + 6, by = "day")) == "Saturday") - 1), "%Y-%m-%d")
  
  
  # Generate a sequence of these end-of-week days
  return(seq.Date(from = first_eow, to = end_date, by = "2 months"))
}

# Generate the specific break points (e.g., all Saturdays)
saturday_breaks <- get_eow_dates(as.Date(min(data_stacked$collection_date), origin=excel_origin), as.Date(max(data_stacked$collection_date), origin=excel_origin))
                                 # day_of_week = "Saturday"))
print(saturday_breaks)

week_ending_dates <- as.Date(c(start), "%Y-%m-%d")

# data_stacked$collection_date <- as.numeric(data_stacked$collection_date)
print(data_stacked$collection_date)

for (date in data_stacked$collection_date) {
  print(date)
  day <- as.Date(date, origin = excel_origin)
  # print(as.Date(as.Date(date, "%Y-%m-%d") + 6, "%Y-%m-%d"))
  week_ending_date_var <- as.Date(day + which(weekdays(seq.Date(day, day + 6, by = "day")) == "Saturday"), "%Y-%m-%d")
  print(week_ending_date_var)
  week_ending_dates <- c(week_ending_dates, week_ending_date_var)
}



# print(week_ending_dates)
data_stacked$week_ending_date <- week_ending_dates[-1]
# print(data_stacked$week_ending_date)
# # data_stacked <- data_stacked %>% fill(week_ending_date, .direction = "down")


# data_stacked <- data_stacked %>% 
#                 filter(week_ending_date == "NA")

# data_stacked$week_ending_date <- as.Date(data_stacked$week_ending_date, "%Y-%m-%d")
print(data_stacked)

# lim <- as.numeric(c(start, end))
# str(lim)
str(saturday_breaks)
# Create the stacked lollipop plot
p <- ggplot(data_stacked, aes(x = week_ending_date,
                              y = position,
                              fill = genotype)) +
  #shape = exposure)) +
  # Add vertical lines (stems) from y=0 to each point
  geom_segment(aes(x = week_ending_date, 
                   xend = week_ending_date,
                   y = 0, 
                   yend = position),
               color = "gray70",
               linewidth = 0.8) +
  # Add points at the top of each stem
  geom_point(size = 5, alpha = 0.8, pch = 21) +
  #scale_shape_manual(values = c("cattle" = 16, 
  #                             "poultry" = 17, 
  #                            "unknown" = 18),
  #                name = "Exposure") +
  scale_fill_manual(values = cols2,
                    #breaks=
                    name = "Genotype") +
  guides(fill = guide_legend(nrow = 1)) +
  #facet_grid(exposure~., ) +
  facet_wrap(~exposure, ncol = 1, scales = "free_x", strip.position = "right") +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  labs(title = "Confirmed H5N1 Feline Cases in North America",
       subtitle = "Nov 2021 - Aug 2025",
       #  x = "Week Ending Date")+
       y = "Weekly Number of Cases", ) +
  theme_minimal(base_size = 10) +
  coord_cartesian(clip = "off") +
  scale_x_date(date_labels = "%b %d, %Y", limits = c(start, end), breaks = saturday_breaks) +
  theme(legend.position = "top", axis.line.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x =element_blank())

# Display the plot
print(p)

# Save the plot
ggsave("C:/Users/maksi/Documents/Statistics/Projects/Avian_Flu_Files/cat_stacked_lollipop.png", 
       plot = p, 
       width = 14, 
       height = 7, 
       dpi = 300)

print("Lollipop plot saved to cat_stacked_lollipop.png")


# Print summary statistics
print("\nSummary by exposure type:")
print(table(data$exposure))
print("\nSummary by genotype:")
print(table(data$genotype))
