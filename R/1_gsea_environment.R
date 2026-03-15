#' @title 交互式构建 GSEA 基因集与超级标签 (Pro 引擎)
#' @description 基于 msigdbr 提供交互式菜单，选择基因集并生成智能语义化标签。支持未知未来子集的动态截断。
#' @param species 物种，默认 "Homo sapiens"
#' @param auto_select 可选。如果不想交互，可传入所选序号向量，例如 c(17, 25)
#' @return 返回包含 TERM2GENE、TERM2NAME、SuperTag 的列表对象
#' @export
build_gsea_pathways_pro <- function(species = "Homo sapiens", auto_select = NULL) {
  # (直接复制您 txt 中 build_gsea_pathways_pro 的原始代码，此函数无需任何改动)
  # ... 由于篇幅，请直接从 txt 拷贝这个函数 ...
}

#' @title 组装计算胶囊 (Pro 引擎)
#' @description 将差异结果(fit)、基因集与表达矩阵完美焊死在一起，实现一次打包，终身复现。
#' @param fit limma 分析得到的 MArrayLM 对象
#' @param pathway_obj build_gsea_pathways_pro() 返回的基因集对象
#' @param expr_data 你的 DGEList 或者标准化后的表达矩阵 (用于画热图，可为 NULL)
#' @export
setup_gsea_env_pro <- function(fit, pathway_obj, expr_data = NULL) {
  if (!inherits(fit, "MArrayLM")) stop("❌ fit 必须是 limma 的对象！")

  c_names <- colnames(fit)
  parsed <- lapply(c_names, function(x) {
    p <- strsplit(x, "\\s*-\\s*")[[1]]
    if(length(p) == 2) c(p[1], p[2]) else c(x, "Unknown")
  })
  parsed_df <- do.call(rbind, parsed)
  contrasts_df <- tibble::tibble(ID = seq_along(c_names), Contrast_Name = c_names, Num = parsed_df[,1], Den = parsed_df[,2])

  env_obj <- list(
    fit = fit,
    contrasts = contrasts_df,
    geneset = list(
      name = pathway_obj$SuperTag,
      term2gene = pathway_obj$TERM2GENE,
      meta_dict = pathway_obj$meta_dict,
      used_collections = pathway_obj$collections_used
    ),
    expr_data = expr_data
  )

  class(env_obj) <- "GseaEnvPro"
  message(sprintf("✅ 胶囊封装完毕！对比组数: %d | Tag: [%s]", nrow(contrasts_df), pathway_obj$SuperTag))
  return(env_obj)
}

#' @title 检查计算胶囊 (Pro 引擎)
#' @param env_pro setup_gsea_env_pro 创建的对象
#' @export
inspect_gsea_env_pro <- function(env_pro) {
  # (直接复制 txt 中的 inspect_gsea_env_pro)
}
