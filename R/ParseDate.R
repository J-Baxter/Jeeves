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

  n <- length(date_str)
  fmt <- rep(NA_character_, n)

  # format detection
  fmt[grepl("^\\d{4}-\\d{2}-\\d{2}$", date_str)] <- "yyyy-mm-dd"
  fmt[grepl("^\\d{2}-\\d{2}-\\d{4}$", date_str)] <- "dd-mm-yyyy"
  fmt[grepl("^\\d{4}-\\d{2}$", date_str)]        <- "yyyy-mm"
  fmt[grepl("^\\d{2}-\\d{4}$", date_str)]        <- "mm-yyyy"
  fmt[grepl("^\\d{4}$", date_str)]               <- "yyyy"

  # initialise
  year  <- month <- day <- rep(NA_integer_, n)

  # yyyy-mm-dd
  i <- fmt == "yyyy-mm-dd"
  if (any(i)) {
    year[i]  <- as.integer(substr(date_str[i], 1, 4))
    month[i] <- as.integer(substr(date_str[i], 6, 7))
    day[i]   <- as.integer(substr(date_str[i], 9, 10))
  }

  # dd-mm-yyyy
  i <- fmt == "dd-mm-yyyy"
  if (any(i)) {
    day[i]   <- as.integer(substr(date_str[i], 1, 2))
    month[i] <- as.integer(substr(date_str[i], 4, 5))
    year[i]  <- as.integer(substr(date_str[i], 7, 10))
  }

  # yyyy-mm
  i <- fmt == "yyyy-mm"
  if (any(i)) {
    year[i]  <- as.integer(substr(date_str[i], 1, 4))
    month[i] <- as.integer(substr(date_str[i], 6, 7))
  }

  # mm-yyyy
  i <- fmt == "mm-yyyy"
  if (any(i)) {
    month[i] <- as.integer(substr(date_str[i], 1, 2))
    year[i]  <- as.integer(substr(date_str[i], 4, 7))
  }

  # yyyy
  i <- fmt == "yyyy"
  if (any(i)) {
    year[i] <- as.integer(date_str[i])
  }

  # build Date (safe defaults)
  date_parsed <- as.Date(
    sprintf(
      "%04d-%02d-%02d",
      year,
      ifelse(is.na(month), 1L, month),
      ifelse(is.na(day),   1L, day)
    )
  )

  out <- data.frame(
    fmt_detected = fmt,
    date_ymd = ifelse(!is.na(day),
                      format(date_parsed, "%Y-%m-%d"),
                      NA_character_),
    date_ym  = ifelse(!is.na(month),
                      format(date_parsed, "%Y-%m"),
                      NA_character_),
    date_y   = ifelse(!is.na(year),
                      format(date_parsed, "%Y"),
                      NA_character_),
    stringsAsFactors = FALSE
  )

  split(out, seq_len(n))
}
