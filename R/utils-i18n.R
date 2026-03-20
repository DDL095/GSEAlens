#' @title GSEAlens 国际化 (i18n) 工具函数
#' @description 管理多语言支持，支持动态语言切换和字典加载。
#' @name utils-i18n
NULL

# 全局语言配置
.GSEALENS_DEFAULT_LANG <- "zh"
.GSEALENS_VALID_LANGS <- c("zh", "en")

#' @title 设置 GSEAlens 语言
#' @description 设置当前会话的语言环境，影响所有 Shiny UI 和消息显示。
#' @param lang 字符串，语言代码。可选 "zh"（中文）或 "en"（英文）。
#' @return 无返回值，设置全局选项。
#' @export
#' @examples
#' \dontrun{
#' set_gsealens_lang("zh")  # 切换到中文
#' set_gsealens_lang("en")  # 切换到英文
#' }
set_gsealens_lang <- function(lang = "zh") {
  if (!lang %in% .GSEALENS_VALID_LANGS) {
    warning(sprintf("不支持的语言代码 '%s'，回退到 '%s'。支持的语言: %s",
                    lang, .GSEALENS_DEFAULT_LANG,
                    paste(.GSEALENS_VALID_LANGS, collapse = ", ")))
    lang <- .GSEALENS_DEFAULT_LANG
  }

  options(GSEALENS.LANG = lang)

  # 同时设置 R 的基础语言（影响系统消息）
  if (lang == "zh") {
    Sys.setenv(LANGUAGE = "zh_CN")
  } else {
    Sys.setenv(LANGUAGE = "en")
  }

  message(sprintf("🌐 GSEAlens 语言已切换为: %s",
                  ifelse(lang == "zh", "中文 (Chinese)", "English")))

  invisible(lang)
}

#' @title 获取当前语言设置
#' @description 内部函数，获取当前设置的语言代码。
#' @return 字符串，当前语言代码。
#' @keywords internal
.get_current_lang <- function() {
  lang <- getOption("GSEALENS.LANG")
  if (is.null(lang)) {
    # 如果未设置，返回默认值并静默设置
    options(GSEALENS.LANG = .GSEALENS_DEFAULT_LANG)
    return(.GSEALENS_DEFAULT_LANG)
  }
  return(lang)
}

#' @title 加载语言字典
#' @description 从 inst/config 目录加载指定语言的 JSON 字典文件。
#' @param lang 语言代码。
#' @return 命名列表，包含所有翻译键值对。
#' @keywords internal
.load_i18n_dict <- function(lang) {
  # 构建文件路径
  file_name <- sprintf("i18n_%s.json", lang)

  # 首先尝试系统文件路径（已安装包）
  file_path <- system.file("config", file_name, package = "GSEAlens")

  # 如果系统路径不存在（开发环境），尝试本地路径
  if (file_path == "" || !file.exists(file_path)) {
    # 尝试从当前工作目录向上查找 inst/config
    possible_paths <- c(
      file.path("inst", "config", file_name),
      file.path("..", "inst", "config", file_name),
      file.path("..", "..", "inst", "config", file_name)
    )

    for (path in possible_paths) {
      if (file.exists(path)) {
        file_path <- path
        break
      }
    }
  }

  # 如果仍未找到，返回空列表
  if (file_path == "" || !file.exists(file_path)) {
    warning(sprintf("找不到语言文件: %s，使用默认键名。", file_name))
    return(list())
  }

  # 读取 JSON 文件
  tryCatch({
    jsonlite::fromJSON(file_path, simplifyVector = TRUE)
  }, error = function(e) {
    warning(sprintf("解析语言文件 %s 失败: %s", file_path, e$message))
    list()
  })
}

#' @title 翻译函数
#' @description 根据键名返回当前语言的翻译文本。支持嵌套键（使用点号分隔）。
#' @param key 字符串，翻译键名。支持嵌套键如 "ui.label_contrast"。
#' @param lang 可选，强制指定语言代码。默认为 NULL，使用当前设置。
#' @param default 可选，如果键不存在时返回的默认值。默认为键名本身。
#' @return 字符串，翻译后的文本。
#' @export
#' @examples
#' \dontrun{
#' .tr("ui.btn_confirm")  # 返回"确认配置"（中文）或"Confirm"（英文）
#' .tr("table.no_data", default = "No data")  # 带默认回退
#' }
.tr <- function(key, lang = NULL, default = NULL) {
  # 确定语言
  if (is.null(lang)) {
    lang <- .get_current_lang()
  }

  # 加载字典（使用缓存避免重复读取文件）
  dict <- .load_i18n_dict_cached(lang)

  # 处理嵌套键（如 "ui.label_contrast"）
  keys <- strsplit(key, "\\.")[[1]]
  result <- dict

  for (k in keys) {
    if (!is.list(result) || is.null(result[[k]])) {
      # 键不存在，返回默认值或键名
      if (!is.null(default)) return(default)
      # 尝试其他语言作为回退
      if (lang != .GSEALENS_DEFAULT_LANG) {
        fallback_result <- .tr(key, lang = .GSEALENS_DEFAULT_LANG, default = key)
        if (fallback_result != key) return(fallback_result)
      }
      return(key)
    }
    result <- result[[k]]
  }

  # 确保返回字符串
  if (is.character(result)) {
    return(result)
  } else {
    return(key)
  }
}

# 字典缓存机制（避免每次翻译都读取文件）
.i18n_cache <- new.env(parent = emptyenv())

#' @title 带缓存的字典加载
#' @keywords internal
.load_i18n_dict_cached <- function(lang) {
  cache_key <- paste0("dict_", lang)

  if (exists(cache_key, envir = .i18n_cache)) {
    return(get(cache_key, envir = .i18n_cache))
  }

  dict <- .load_i18n_dict(lang)
  assign(cache_key, dict, envir = .i18n_cache)
  return(dict)
}

#' @title 清除翻译缓存
#' @description 清除内部字典缓存，用于开发调试或热重载语言文件。
#' @export
clear_i18n_cache <- function() {
  rm(list = ls(envir = .i18n_cache), envir = .i18n_cache)
  message("🌐 语言缓存已清除")
  invisible(NULL)
}

#' @title 获取所有可用翻译键
#' @description 用于开发调试，列出当前语言的所有可用键。
#' @param lang 语言代码，默认当前语言。
#' @return 字符向量，所有键名（展平格式，使用点号分隔）。
#' @export
list_i18n_keys <- function(lang = NULL) {
  if (is.null(lang)) lang <- .get_current_lang()
  dict <- .load_i18n_dict(lang)

  # 递归展平列表
  flatten_keys <- function(x, prefix = "") {
    if (!is.list(x)) return(prefix)

    keys <- names(x)
    if (is.null(keys)) return(prefix)

    result <- character(0)
    for (k in keys) {
      new_prefix <- if (prefix == "") k else paste(prefix, k, sep = ".")
      if (is.list(x[[k]]) && length(x[[k]]) > 0 && !is.null(names(x[[k]]))) {
        result <- c(result, flatten_keys(x[[k]], new_prefix))
      } else {
        result <- c(result, new_prefix)
      }
    }
    return(result)
  }

  sort(flatten_keys(dict))
}
