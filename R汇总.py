from pathlib import Path
from datetime import datetime


def read_file_with_fallback(file_path: Path) -> str:
    """
    尝试用多种常见编码读取文件内容。
    如果都失败，则使用 utf-8 并替换非法字符。
    """
    encodings = ["utf-8", "utf-8-sig", "gb18030", "gbk", "latin-1"]

    for enc in encodings:
        try:
            return file_path.read_text(encoding=enc)
        except UnicodeDecodeError:
            continue
        except Exception as e:
            return f"[读取失败] {file_path}\n错误信息: {e}\n"

    # 最后兜底
    try:
        return file_path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return f"[读取失败] {file_path}\n错误信息: {e}\n"


def main():
    # 脚本所在目录
    script_dir = Path(__file__).resolve().parent

    # 递归查找所有 .R / .r 文件
    r_files = sorted(
        [p for p in script_dir.rglob("*") if p.is_file() and p.suffix.lower() == ".r"],
        key=lambda x: str(x.relative_to(script_dir)).lower()
    )

    # 时间戳文件名
    timestamp = datetime.now().strftime("%Y-%m-%d-%H%M%S")
    output_file = script_dir / f"{timestamp} GSEAlens全量R代码.txt"

    with output_file.open("w", encoding="utf-8") as out:
        out.write(f"GSEAlens 全量 R 代码汇总\n")
        out.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        out.write(f"扫描目录: {script_dir}\n")
        out.write(f"R 文件数量: {len(r_files)}\n")
        out.write("=" * 100 + "\n\n")

        if not r_files:
            out.write("未找到任何 .R 文件。\n")
        else:
            for i, file_path in enumerate(r_files, start=1):
                relative_path = file_path.relative_to(script_dir)
                content = read_file_with_fallback(file_path)

                out.write("=" * 100 + "\n")
                out.write(f"[文件 {i}/{len(r_files)}] {relative_path}\n")
                out.write("=" * 100 + "\n")
                out.write(content)

                if not content.endswith("\n"):
                    out.write("\n")
                out.write("\n\n")

    print(f"汇总完成：{output_file}")
    print(f"共处理 {len(r_files)} 个 R 文件。")


if __name__ == "__main__":
    main()
