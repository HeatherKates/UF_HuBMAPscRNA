mod_labels <- function(object, values, select = c("first", "last", "common", "all"), simplify = TRUE, ...) {
  select <- select[1L]
  select <- match.arg(arg = select)
  values <- intersect(values, rownames(object))
  
  if (length(values) == 0) {
    return(character(0))
  }
  
  # Precompute logical matrix and column names
  cmap_data <- as.matrix(object)
  colnames_object <- colnames(object)
  rownames_object <- rownames(object)
  
  # Get row indices for values
  row_indices <- match(values, rownames_object)
  
  # Initialize the list to store results
  obs <- vector("list", length(values))
  names(obs) <- values
  
  # Direct indexing to replace sapply
  for (i in seq_along(row_indices)) {
    row_idx <- row_indices[i]
    vals <- colnames_object[cmap_data[row_idx, , drop = FALSE]]
    if (length(vals) > 0) {
      obs[[i]] <- vals
    }
  }
  
  obs <- Filter(length, obs)
  
  obs <- switch(select, 
                first = lapply(obs, `[[`, 1L), 
                last = lapply(obs, function(x) x[[length(x)]]), 
                common = {
                  counts <- table(unlist(obs))
                  tmp <- obs
                  obs <- vector("character", length(tmp))
                  names(obs) <- names(tmp)
                  for (i in seq_along(obs)) {
                    obs[i] <- names(which.max(counts[names(counts) %in% tmp[[i]]]))
                  }
                  obs
                }, 
                obs)
  
  if (isTRUE(simplify)) {
    tmp <- obs
    obs <- unlist(tmp)
    names(obs) <- make.unique(rep(names(tmp), times = lengths(tmp)))
  }
  
  return(obs)
}
