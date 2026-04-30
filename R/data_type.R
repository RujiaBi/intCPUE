.detect_data_type_intCPUE <- function(data_utm, data_type = c("auto", "delta", "positive", "mixed")) {
  data_type <- match.arg(data_type)
  .check_required_cols(data_utm, c("encounter", "flagid"))
  if (anyNA(data_utm$encounter) || anyNA(data_utm$flagid)) {
    stop("`encounter` and `flagid` must not contain NA.", call. = FALSE)
  }
  if (!all(data_utm$encounter %in% c(0L, 1L))) {
    stop("`encounter` must contain only 0/1.", call. = FALSE)
  }

  by_flag <- split(as.integer(data_utm$encounter), data_utm$flagid)
  has_zero <- vapply(by_flag, function(x) any(x == 0L), logical(1))
  has_pos <- vapply(by_flag, function(x) any(x == 1L), logical(1))

  if (any(!has_pos)) {
    bad <- names(has_pos)[!has_pos]
    stop(
      "Each flag/type must have at least one positive observation. Flags with only zero catch: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  flag_type <- ifelse(has_zero & has_pos, "delta", "positive")
  detected <- if (all(flag_type == "delta")) {
    "delta"
  } else if (all(flag_type == "positive")) {
    "positive"
  } else {
    "mixed"
  }

  final <- if (identical(data_type, "auto")) detected else data_type
  if (!identical(data_type, "auto") && !identical(final, detected)) {
    stop(
      "`data_type='", data_type, "'` is inconsistent with the data; auto-detected `",
      detected, "`. Use `data_type='auto'` or check encounter by flag/type.",
      call. = FALSE
    )
  }

  flag_ids <- suppressWarnings(as.integer(names(flag_type)))
  if (anyNA(flag_ids)) {
    flag_ids <- seq_along(flag_type) - 1L
  }
  ord <- order(flag_ids)
  flag_ids <- flag_ids[ord]
  flag_type <- unname(flag_type[ord])

  if (identical(final, "mixed")) {
    ref_pos <- match(0L, flag_ids)
    if (is.na(ref_pos)) {
      stop("`data_type='mixed'` requires reference flag `flagid = 0`.", call. = FALSE)
    }
    if (!identical(flag_type[ref_pos], "delta")) {
      stop(
        "`data_type='mixed'` requires the reference flag (`flagid = 0`) to be delta ",
        "(i.e., to include both zero and positive observations).",
        call. = FALSE
      )
    }
  }

  list(
    requested = data_type,
    data_type = final,
    detected = detected,
    flag_id = flag_ids,
    flag_type = flag_type,
    flag_uses_component1 = flag_type == "delta"
  )
}

.validate_delta_link_intCPUE <- function(delta_link = c("Poisson", "traditional"), data_type) {
  delta_link <- match.arg(delta_link)
  if (identical(data_type, "positive")) {
    return("traditional")
  }
  if (identical(data_type, "mixed") && identical(delta_link, "Poisson")) {
    stop(
      "`delta_link='Poisson'` is only valid when `data_type='delta'`. ",
      "`data_type='mixed'` requires `delta_link='traditional'` because ",
      "positive-only flags have encounter probability fixed at 1.",
      call. = FALSE
    )
  }
  delta_link
}
