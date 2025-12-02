# Load required libraries
library(readxl)
library(ggplot2)
library(dplyr)

# Read the data
data <- read_excel("/mnt/user-data/uploads/all-2344b-seqs-gisaid-metadata.xlsx")


# Preview the data
print("Data summary:")
print(head(data))
print(paste("Total rows:", nrow(data)))

# Remove rows with NA in week ending date
data <- data %>% 
  filter(!is.na(`week ending date`))

# For stacked lollipops, we need to count occurrences per week/genotype/exposure combination
# and create a position variable for stacking
data_stacked <- data %>%
  arrange(`week ending date`, exposure, genotype) %>%
  group_by(`week ending date`, exposure) %>%
  mutate(position = row_number()) %>%
  ungroup()


#cols<- c("B3.13"= "orange", "B3.13 (inferred)"= "#FED8B1", "D1.1"= "#7F00FF", "D1.1 (inferred)"= "#D6B4FC",
 #        "B3.2"= "#895129", "D1.3"= "steelblue", "unknown"= "grey")

cols2<- c("B3.13"= "orange", "B3.13 (inferred)"= "#FED8B1", "D1.1"= "#7F00FF", "D1.1 (inferred)"= "#D6B4FC",
          "D1.3"= "steelblue", "unknown"= "grey")

start<- as.Date("2024-03-01")
end<- as.Date("2025-03-23")

data_stacked$genotype<- factor(data_stacked$genotype, levels = c(""))



# Function to find the next desired end-of-week date (e.g., the first Saturday after a start date)
get_eow_dates <- function(start_date, end_date, day_of_week = "Saturday") {
  # Convert start and end dates to Date objects if they aren't already
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  
  # Find the first end-of-week day (e.g., Saturday) on or after the start date
  first_eow <- start_date + (which(weekdays(seq(start_date, start_date + 6, by = "day")) == day_of_week) - 1)
  
  # Generate a sequence of these end-of-week days
  seq.Date(from = first_eow, to = end_date, by = "2 week")
}

# Generate the specific break points (e.g., all Saturdays)
saturday_breaks <- get_eow_dates(min(data_stacked$`week ending date`), max(data_stacked$`week ending date`),
                                 day_of_week = "Saturday")

# Create the stacked lollipop plot
p <- ggplot(data_stacked, aes(x = `week ending date`, 
                              y = position,
                              fill = genotype))+
  #shape = exposure)) +
  # Add vertical lines (stems) from y=0 to each point
  geom_segment(aes(x = `week ending date`, 
                   xend = `week ending date`,
                   y = 0, 
                   yend = position),
               color = "gray70",
               linewidth = 0.8) +
  # Add points at the top of each stem
  geom_point(size = 5, alpha = 0.8, pch= 21) +
  #scale_shape_manual(values = c("cattle" = 16, 
  #                             "poultry" = 17, 
  #                            "unknown" = 18),
  #                name = "Exposure") +
  scale_fill_manual(values = cols2,
                    #breaks=
                    name = "Genotype") +
  guides(fill = guide_legend(nrow = 1))+
  #facet_grid(exposure~., ) +
  facet_wrap(~exposure, ncol = 1, scales = "free_x", strip.position = "right")+
  geom_hline(yintercept=0, color = "black", size=1) +
  labs(title = "Confirmed H5N1 Cases in North America",
       subtitle = "Mar 2024 - Mar 2025",
       #  x = "Week Ending Date")+
       y = "Weekly Number of Cases", ) +
  theme_minimal(base_size = 10) +
  coord_cartesian(clip = "off")+
  scale_x_date(date_labels = "%b %d, %Y", limits=c(start,end), breaks = saturday_breaks)+
  theme(legend.position = "top", axis.line.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x =element_blank())

# Display the plot
print(p)

# Save the plot
ggsave("/mnt/user-data/outputs/h5n1_stacked_lollipop.png", 
       plot = p, 
       width = 14, 
       height = 7, 
       dpi = 300)

print("Lollipop plot saved to /mnt/user-data/outputs/h5n1_stacked_lollipop.png")


# Print summary statistics
print("\nSummary by exposure type:")
print(table(data$exposure))
print("\nSummary by genotype:")
print(table(data$genotype))