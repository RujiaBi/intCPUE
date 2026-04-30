.template_id_intCPUE <- function(q_diffs_spatial = c("on", "off"),
                                 pop_spatial = c("on", "off"),
                                 pop_spatiotemporal = c("on", "off")) {
  q_diffs_spatial <- match.arg(q_diffs_spatial)
  pop_spatial <- match.arg(pop_spatial)
  pop_spatiotemporal <- match.arg(pop_spatiotemporal)
  omega_bit <- as.integer(identical(pop_spatial, "on"))
  epsilon_bit <- as.integer(identical(pop_spatiotemporal, "on"))
  flag_s_bit <- as.integer(identical(q_diffs_spatial, "on"))
  paste0("intCPUE_w", omega_bit, epsilon_bit, flag_s_bit)
}

.ensure_template_dll_intCPUE <- function(q_diffs_spatial = c("on", "off"),
                                         pop_spatial = c("on", "off"),
                                         pop_spatiotemporal = c("on", "off")) {
  template_id <- .template_id_intCPUE(
    q_diffs_spatial = q_diffs_spatial,
    pop_spatial = pop_spatial,
    pop_spatiotemporal = pop_spatiotemporal
  )

  dll_dir <- system.file("tmb", package = "intCPUE")
  if (!nzchar(dll_dir)) {
    stop("Could not locate installed intCPUE template DLLs in inst/tmb.", call. = FALSE)
  }

  dll_path <- TMB::dynlib(file.path(dll_dir, template_id))
  if (!file.exists(dll_path)) {
    stop(
      "Missing precompiled template DLL for ", template_id,
      ". Reinstall intCPUE so installation can build the template set.",
      call. = FALSE
    )
  }

  loaded_names <- names(getLoadedDLLs())
  if (!(template_id %in% loaded_names)) {
    dyn.load(dll_path)
  }

  template_id
}
