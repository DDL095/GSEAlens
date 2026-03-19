用于多组比较的GSEA数据的探索与发现。名称中的lens代表这个包如同放大透镜一样，探索相关的gsea通路。
通过在网页内展示通路介绍与通路说明，加入便于AI导出通路富集情况，让GSEA富集分析更加方便与便捷，并且还可以根据所需信息，导入经过翻译后的通路注释与摘要说明。
这个R包通过打包流程，规范输入内容，能够简化GSEA富集分析的查看与探索过程。
由于多组比较在limma-voom流程中使用无截距组别设置会更加方便，本项目的代码基于无截距组别设置进行结果分析与处理。

基于参考文献*****，GSEAlens 将接受deseq2::result中的wald值(result$stat)与limma::fit的t值(fit$t)作为fgsea的输入并进行排序，对应的比对将基于输入数据中包含的分组信息自动生成。因此，用户需要指定的目标分类因子（例如 condition 或 group），deseq2流程允许模型中包含额外的加性协变量（例如 batch、subject、sex）用于校正。在limma-voom流程中，需要指定分组为无截距方式(~0+group);在deseq2当中支持~ group，~ batch + group，~ subject + group，但最终均只支持针对group的分组分析。GSEAlens 暂不考虑支持交互项、连续变量效应或复杂自定义 contrast 的批量 GSEA 比对结果的生成与探索。
