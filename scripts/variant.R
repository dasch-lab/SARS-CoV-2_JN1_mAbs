#!/usr/bin/Rscript

usePackage <- function(package){
  if(eval(parse(text=paste("require(",package,")")))) {
    return (TRUE)
  }

  install.packages(package, repos='http://cran.us.r-project.org')
  return (eval(parse(text=paste("require(",package,")"))))
}

usePackage.bio <- function(package){

  library.path <- .libPaths()[1];

  if(!is.element(package, installed.packages(lib.loc=library.path)[,1]))
    BiocManager::install(package, lib=library.path, update=FALSE)

  return (eval(parse(text=paste("require(",package,")"))))
}

# For command line parsing
usePackage("argparse");

# create parser object
parser <- ArgumentParser()
parser$add_argument("--frequency", type="character", default=NULL, help="input GISAID frequency tsv file");
parser$add_argument("--functionality", type="character", default=NULL, help="input MADLAB functionality tsv file");
parser$add_argument("-o", "--output", default=NULL, type="character", help="output pdf file");

opt <- parser$parse_args();

# Check options
if (is.null(opt$frequency) || is.null(opt$output) || is.null(opt$functionality)){
  stop("Missing argument(s)!", call.=FALSE)
}

# Get output path and directory
output_directory <- dirname(opt$output);
output_name <- tools::file_path_sans_ext(basename(opt$output));
output_extension <- tools::file_ext(opt$output);

# Install packages
usePackage("ggplot2");
usePackage("tidyr");
usePackage('gridExtra');
usePackage("dplyr");
usePackage("cowplot");  # For plot_grid
usePackage("RColorBrewer");
usePackage("ggnewscale");
usePackage("stringr");
usePackage("ggrepel");
usePackage("ggnewscale");

# Set maximum number of rows to be displayed
options(dplyr.print_max = 1E9);

# Function to convert date format
convert_date <- function(date_list) {

  result <- c();
  for(date_string in date_list){
    # Convert string to Date object
    date_formatted <- paste(substr(date_string, 1, 4), substr(date_string, 5, nchar(date_string)), sep = "-");
    date <- as.Date(paste(date_formatted, "01", sep="-"));

    # Format the date to "Month Year"
    formatted_date <- format(date, "%b %Y");
    result <- c(result, formatted_date);
  }

  return(result);
}

# Load dataset
data.frequency <- read.csv(opt$frequency, header = TRUE, sep="\t", stringsAsFactors=FALSE) %>%
  arrange(date) %>%
  mutate(id = match(date, unique(date))-1) %>%
  mutate(date = convert_date(date)) %>%
  mutate(freq = 100*freq);

# Load functionality
data.functionality <- read.csv(opt$functionality, header=TRUE, sep='\t', stringsAsFactors=FALSE) %>%
  dplyr::mutate(value=ifelse(value < 100000, 1, 0)) %>%
  dplyr::group_by(germline, variant) %>%
  dplyr::summarise(activity=round(100*sum(value)/n(),2));

# Create a palette
# https://coolors.co/ffe873-d7b40f-af8030-bf5656
variant_list <- unique(data.frequency %>% pull(variant));
variant_color <- setNames(gray.colors(length(variant_list), start = 0.5, end = 0.9), variant_list)
selected_germlines <- c('IGHV3-53;IGHJ6-1'='#FFE873', 'IGHV3-66;IGHJ4-1'='#D7B40F', 'IGHV3-66;IGHJ6-1'='#AF8030', 'IGHV5-10-1;IGHJ4-1'='#BF5656');
selected_variant <- c('XBB.1.5'='#BEDACD','EG.5.1.1'='#298072', 'BA.2.86'='#E8A2FF', 'JN.1'='#D14AFD');
for (key in names(selected_variant)) {
  variant_color[key] <- selected_variant[key];  # Set value by key to another vector
}

cat('Germline functionality:\n')
print(data.functionality %>% filter(germline %in% names(selected_germlines)) %>% filter(variant %in% names(selected_variant)))

# Convert to factors
data.frequency$variant <- factor(variant_list, levels = rev(c(names(selected_variant), setdiff(unique(variant_list), names(selected_variant)))))
data.frequency$date <- factor(data.frequency$date , levels=sort(unique(data.frequency$date)) );

# Calculate adjusted frequency for each germline
data.germline <- NULL;
germline_list <- data.functionality %>% distinct(germline) %>% pull();
id_list <- data.frequency %>% distinct(id) %>% arrange(id) %>% pull();
for (key in names(selected_germlines)){

  data.germline.single <- data.frame(
    id = integer(), 
    germline = character(), 
    activity = double()
  );

  data.variant <- data.functionality %>%
    filter(germline == key) %>%
    filter(variant %in% names(selected_variant));
  if (nrow(data.variant) == 0){
    next;
  }
  for (date_id in id_list) {
    date_value <- 0;
    for (i in 1:nrow(data.variant)){
      row <- data.variant[i, ];
      variant_value <- data.frequency %>%
        filter(id == date_id) %>%
        filter(variant == row$variant) %>%
        pull(freq);
      if (length(variant_value) == 0){
        next;
      }

      # Sum functionality over variants
      date_value <- date_value + (row$activity * variant_value/100);
    }

    # Append row
    data.germline.single <- rbind(data.germline.single, 
      data.frame(
        id=date_id, 
        germline=row$germline, 
        activity=date_value
      )
    );
  
  }

  # Store data
  if (length(data.germline) == 0){
    data.germline <- data.germline.single;
  }else{
    data.germline <- rbind(data.germline,data.germline.single); 
  }
}

# Calculate log scale
data.germline <- data.germline %>%
  mutate(activity_log = log10(activity+1));

# Convert to wide for supplementary
data.germline.wide <- data.germline %>%
  select(id, activity, germline) %>%
  pivot_wider(names_from = id, values_from = activity);
write.table(data.germline.wide, file = file.path(output_directory, "variant_germline.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Ratio of largest NewCases value to largest activity_log10 value
# Also multiplied by 1.5 to make the grey area plot extend above the bars
ratio = 1.0 + 87/max(data.germline$activity_log);

data.border <- NULL;
for (date_id in id_list){
  cumulative <- 0;
  for (key in names(selected_variant)) {
    variant_value <- data.frequency %>%
        filter(id == date_id) %>%
        filter(variant == key) %>%
        pull(freq);
    cumulative <- cumulative + variant_value;
    data.border.single <- data.frame(
        id=date_id, 
        variant=key, 
        freq=cumulative
      );
    if (length(data.border) == 0){
      data.border <- data.border.single;
    }else{
      data.border <- rbind(data.border, data.border.single); 
    }
  }
}

# Plot results
p.variant <- ggplot() +
  geom_area(data.frequency, mapping=aes(x = id, y = freq, fill = variant)) + 
  geom_line(data=data.border, mapping=aes(x=id, y=freq, group=variant), colour='black', linewidth=0.2) +
  scale_fill_manual(values=variant_color, breaks = names(selected_variant),guide = "none") + 
  new_scale_fill() + 
  geom_line(data = data.germline, mapping=aes(x = id, y = activity_log*ratio, colour=germline)) +
  geom_point(data = data.germline %>% filter(!id %in% c(0, 12)), mapping=aes(x = id, y = activity_log*ratio, fill=germline),  colour="black",pch=21, size=2) + 
  scale_y_continuous(
    expand=c(0,0), 
    name="Variant frequency (%)",
    sec.axis = sec_axis(~ ./ratio, breaks=c(0:5, 0:5),labels=function(x) round(10^x,2), 
    name='Gemline-descendent mAbs activity (%)'
  )) +
  scale_x_continuous(breaks = unique(data.frequency$id), labels = unique(data.frequency$date), expand=c(0,0),guide = guide_axis(angle = 90)) +
  scale_fill_manual(values=selected_germlines) + 
  scale_colour_manual(values=selected_germlines) + 
  labs(y = "Frequency (%)")  +
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=8, hjust=1,margin = margin(t=0, r=10, b=0, l=0)),
    axis.ticks.x = element_line(color='black'),
    axis.ticks.y = element_line(color='black'),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
  );

# Get legend
p.variant.legend <- get_legend(p.variant + theme(
  legend.title=element_blank(),
  legend.direction = "horizontal",
  legend.justification="center",
  legend.box.just = "bottom", 
  legend.text=element_text(size=6)
));

# Remove the legend before plotting
p.variant <- p.variant + theme(legend.position='none');
  
grid <- plot_grid(
  p.variant,
  p.variant.legend,
  nrow=2,
  ncol=1,
  rel_heights=c(1.0, 0.2)
) + theme(
  plot.margin = margin(t=0.3, r=0.3, b=0.3, l=0.3, unit="cm"),
);

# Save the final plot
ggsave(filename=opt$output, plot=grid, units="mm", width=170, height=120, dpi=270);
ggsave(filename=file.path(output_directory, paste0(output_name, '.',"png")), plot=grid, units="mm", width=170, height=120, dpi=600);

# Often the most useful way to look at many warnings:
summary(warnings())
