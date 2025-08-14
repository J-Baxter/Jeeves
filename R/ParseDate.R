#' ParseDate
#'
#' Format a date string with variable precision/delimiters. Output a
#'
#' @param date_str date string. May take yyyy-mm-dd, dd-mm-yyyy, yyyy-mm, mm-yyyy
#' or yyyy formats. Caution advised if using mm-dd-yyyy.
#'
#' @returns A list-column of formats including ymd, ym, y, decimal and precision.
#' @export
#'
#' @examples df <- tibble(collection_date = c("2025-08-14",
#'  "14-08-2025", "2025-08", "08-2025", "2025"))
#'
#' df_nested <- df %>%
#'  mutate(parsed_info = ParseDate(collection_date)) %>%
#'  unnest(parsed_info)
#'
ParseDate <- function(date_str) {

  # Detect format with grepl (vectorised)
  fmt <- case_when(
    str_detect(date_str, "^\\d{4}-\\d{2}-\\d{2}$") ~ "yyyy-mm-dd",
    str_detect(date_str, "^\\d{2}-\\d{2}-\\d{4}$") ~ "dd-mm-yyyy",
    str_detect(date_str, "^\\d{4}-\\d{2}$")        ~ "yyyy-mm",
    str_detect(date_str, "^\\d{2}-\\d{4}$")        ~ "mm-yyyy",
    str_detect(date_str, "^\\d{4}$")               ~ "yyyy",
    TRUE                                           ~ NA_character_
  )

  # Vectorised parse
  date_parsed <- as.Date(rep(NA_character_, length(date_str)))
  date_parsed[fmt == "yyyy-mm-dd"] <- ymd(date_str[fmt == "yyyy-mm-dd"])
  date_parsed[fmt == "dd-mm-yyyy"] <- dmy(date_str[fmt == "dd-mm-yyyy"])
  date_parsed[fmt == "mm-dd-yyyy"] <- mdy(date_str[fmt == "mm-dd-yyyy"])
  date_parsed[fmt == "yyyy-mm"]    <- ym(date_str[fmt == "yyyy-mm"])
  date_parsed[fmt == "mm-yyyy"]    <- my(date_str[fmt == "mm-yyyy"])
  date_parsed[fmt == "yyyy"]       <- as.Date(as.POSIXlt(date_str[fmt == "yyyy"], format = "%Y"))

  # Return list-column of tibbles (one per row)
  tibble(
    fmt_detected = fmt,
    date_parsed  = date_parsed,
    date_ymd     = ifelse(fmt %in% c("yyyy-mm-dd", "dd-mm-yyyy", "mm-dd-yyyy"),
                          format(date_parsed, "%Y-%m-%d"), NA_character_),
    date_ym      = ifelse(fmt %in% c("yyyy-mm-dd", "dd-mm-yyyy", "mm-dd-yyyy", "yyyy-mm", "mm-yyyy"),
                          format(date_parsed, "%Y-%m"), NA_character_),
    date_y       = ifelse(!is.na(date_parsed), format(date_parsed, "%Y"), NA_character_)
  ) %>%
    split(1:nrow(.)) %>%
    unname()
}
