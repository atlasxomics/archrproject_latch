library("ggplot2")

# Color defaults --------------------------------------------------------------
c_list <- list(
  c(
    "grey85", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84",
    "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"
  ),
  c(
    "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
    "#FEE090", "#FDAE61", "#F46D43", "#D73027"
  )[c(1, 1:9, 9)],
  c(
    "#FDE725", "#AADC32", "#5DC863", "#27AD81", "#21908C",
    "#2C728E", "#3B528B", "#472D7B", "#440154"
  ),
  c(
    "#E6E7E8", "#3A97FF", "#8816A7", "black"
  ),
  c(
    "#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
    "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"
  )
)
names(c_list) <- c(
  "White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple", "comet", "blueYellow"
)

# Panel sizes -----------------------------------------------------------------
p_list <- c("400px", "600px", "800px")
names(p_list) <- c("Small", "Medium", "Large")

p_list2 <- c("500px", "700px", "900px")
names(p_list2) <- c("Small", "Medium", "Large")

p_list3 <- c("600px", "800px", "1000px")
names(p_list3) <- c("Small", "Medium", "Large")

s_list <- c(18, 24, 30)
names(s_list) <- c("Small", "Medium", "Large")

lList <- c(5, 6, 7)
names(lList) <- c("Small", "Medium", "Large")


sctheme <- function(base_size = 24, xy_val = TRUE, x_ang = 0, x_jus_h = 0.5) {
  oup_theme <- theme(
    text = element_text(size = base_size, family = "Helvetica"),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black", size = base_size / 20),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = base_size),
    axis.text.x = element_text(angle = x_ang, hjust = x_jus_h),
    legend.position = "bottom",
    legend.key = element_rect(colour = NA, fill = NA),
  )
  if (!xy_val) {
    oup_theme <- oup_theme + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank()
    )
  }
  return(oup_theme)
}
