#' @title 初始化 GSEA 环境对象
#' @description 统一入口，支持 limma-voom (无截距设计) 和 DESeq2 对象。
#' 自动提取对比组、差异分析结果和表达矩阵，生成标准化的 GseaEnv 对象。
#' @param fit 差异分析对象。支持 \code{MArrayLM} (limma) 或 \code{DESeqDataSet} (DESeq2)。
#' @param pathway_obj 基因集对象。通常由 \code{build_gsea_pathways()} 生成。
#' @param expr_data 可选。表达矩阵或 DGEList。若为 NULL，将尝试从 fit 对象中提取 (仅部分支持)。
#' @param target_factor 仅用于 DESeq2。指定感兴趣的目标因子。若为 NULL，自动推断为设计公式最后一项。
#' @return 返回一个类为 \code{GseaEnv} 的标准化对象。
#' @export
#' @examples
#' \dontrun{
#' # Limma 流程
#' design <- model.matrix(~ 0 + group, data = samples)
#' fit <- lmFit(expr, design) %>% contrasts.fit(cont.matrix) %>% eBayes()
#' env <- setup_gsea_env(fit, pathway_obj)
#'
#' # DESeq2 流程
#' dds <- DESeqDataSetFromMatrix(countData, colData, design = ~ batch + group)
#' dds <- DESeq(dds)
#' env <- setup_gsea_env(dds, pathway_obj, target_factor = "group")
#' }
setup_gsea_env <- function(fit, pathway_obj, expr_data = NULL, target_factor = NULL) {

  message("\n🚀 GSEAlens 引擎启动...")

  # 1. 校验基因集对象
  if (is.null(pathway_obj) || is.null(pathway_obj$TERM2GENE)) {
    stop("❌ pathway_obj 无效！请使用 build_gsea_pathways() 生成有效的基因集对象。")
  }

  # 2. 识别后端类型并分发提取
  backend_data <- NULL
  backend_info <- NULL

  if (inherits(fit, "MArrayLM")) {
    message("🔍 检测到输入类型: [Limma-Voom]")

    # 调用 Limma 提取器
    backend_data <- .extract_limma_data(fit, expr_data)

    backend_info <- list(
      backend = "limma_voom",
      input_class = "MArrayLM",
      design_formula = NULL, # Limma fit 对象通常不直接保留公式
      target_factor = NA
    )

  } else if (inherits(fit, "DESeqDataSet")) {
    message("🔍 检测到输入类型: [DESeq2]")

    # 调用 DESeq2 提取器
    backend_data <- .extract_deseq2_data(fit, target_factor)

    # 记录实际使用的 target_factor
    actual_factor <- if(is.null(target_factor)) {
      utils::tail(attr(terms(DESeq2::design(fit)), "term.labels"), 1)
    } else {
      target_factor
    }

    backend_info <- list(
      backend = "deseq2",
      input_class = "DESeqDataSet",
      design_formula = DESeq2::design(fit),
      target_factor = actual_factor
    )

  } else {
    stop("❌ 不支持的输入类型！请传入 MArrayLM (limma) 或 DESeqDataSet (DESeq2) 对象。")
  }

  # 3. 组装基因集信息
  geneset_info <- list(
    name = pathway_obj$SuperTag,
    term2gene = pathway_obj$TERM2GENE,
    meta_dict = pathway_obj$meta_dict,
    used_collections = pathway_obj$collections_used
  )

  # 4. 构建最终对象
  env_obj <- create_gsea_env(
    backend_info = backend_info,
    contrast_registry = backend_data$contrast_registry,
    de_store = backend_data$de_store,
    expr_bundle = backend_data$expr_bundle,
    geneset = geneset_info,
    raw_obj = list(fit = fit) # 保留原始对象引用
  )

  # 5. 最终校验
  .check_gsea_env(env_obj)

  message("✅ GseaEnv 对象构建成功！")
  message(sprintf("   📦 包含 %d 个对比组", nrow(backend_data$contrast_registry)))
  message(sprintf("   🧬 包含 %d 条通路", nrow(pathway_obj$meta_dict)))

  return(env_obj)
}


# 查看函数


#' @title 查看 GSEA 环境对象概览
#' @description 在控制台打印 GseaEnv 对象的详细信息，包括后端类型、对比组列表和基因集状态。
#' @param env_obj GseaEnv 对象
#' @export
inspect_gsea_env <- function(env_obj) {
  if (!inherits(env_obj, "GseaEnv")) stop("输入对象不是 GseaEnv 类。")

  # 解析信息
  bi <- env_obj$backend_info
  reg <- env_obj$contrast_registry
  gs <- env_obj$geneset
  expr <- env_obj$expr_bundle

  # 控制台输出
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("       🧬 GSEAlens Environment Summary\n")
  cat(rep("=", 60), "\n\n", sep = "")

  # 1. 后端信息
  cat("⚙️  [1] Backend Information\n")
  cat(sprintf("   • Type        : %s\n", bi$backend))
  cat(sprintf("   • Input Class : %s\n", bi$input_class))
  if (bi$backend == "deseq2") {
    cat(sprintf("   • Target Factor: %s\n", bi$target_factor))
  }
  cat("\n")

  # 2. 对比组注册表
  cat("⚖️  [2] Contrast Registry (", nrow(reg), " comparisons)\n", sep = "")
  if (nrow(reg) > 0) {
    # 打印前几个，防止刷屏
    print_head <- utils::head(reg[, c("contrast_id", "left_group", "right_group")], 5)
    for (i in 1:nrow(print_head)) {
      cat(sprintf("   [%d] %s  ( %s vs %s )\n",
                  i, print_head$contrast_id[i],
                  print_head$left_group[i], print_head$right_group[i]))
    }
    if (nrow(reg) > 5) cat("   ... (省略剩余)\n")
  } else {
    cat("   ⚠️ 未检测到对比组！\n")
  }
  cat("\n")

  # 3. 基因集信息
  cat("🧫 [3] Gene Set Database\n")
  cat(sprintf("   • Name        : %s\n", gs$name))
  cat(sprintf("   • Pathways    : %d\n", nrow(gs$meta_dict)))
  cat("\n")

  # 4. 表达数据状态
  cat("📊 [4] Expression Data Status\n")
  if (!is.null(expr$raw_counts)) {
    cat(sprintf("   • Raw Counts  : %d genes x %d samples\n",
                nrow(expr$raw_counts), ncol(expr$raw_counts)))
    cat(sprintf("   • Display Mat : %s\n", ifelse(!is.null(expr$display_expr), "Ready", "Missing")))
  } else {
    cat("   ⚠️ 未包含表达矩阵 (后续热图功能将不可用)\n")
  }
  cat("\n")

  # 5. 下一步指引
  cat("🚀 [5] Next Step\n")
  cat("   请运行计算引擎：\n")
  cat("   > res <- batch_calc_gsea(env_obj)\n")
  cat(rep("=", 60), "\n\n", sep = "")

  invisible(env_obj)
}

#' @export
print.GseaEnv <- function(x, ...) {
  inspect_gsea_env(x)
}
