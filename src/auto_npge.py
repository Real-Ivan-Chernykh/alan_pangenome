import argparse
import subprocess
from datetime import datetime
from pathlib import Path

# Список команд, которые нужно выполнить по порядку
COMMANDS = [
    ["npge", "Prepare"],
    ["npge", "Examine"],
    ["npge", "MakePangenome"],
    ["npge", "CheckPangenome"],
    ["npge", "PostProcessing"],
]

LOG_NAME = "npge_errors.log"


def run_commands(target_dir: str) -> None:
    directory = Path(target_dir).expanduser().resolve()

    # Проверяем, что директория существует. Если нет — выбрасываем ошибку.
    if not directory.is_dir():
        raise FileNotFoundError(f"Директория не найдена: {directory}")
    print(f"Рабочая директория: {directory}")

    log_path = directory / LOG_NAME

    for cmd in COMMANDS:
        print(f"\n=== Выполняется: {' '.join(cmd)} ===")
        result = subprocess.run(
            cmd,
            cwd=directory,   # запускаем команду именно из этой директории
            text=True,
            capture_output=True,   # перехватываем stdout и stderr
        )

        print(
                f"\nКоманда '{' '.join(cmd)}' завершилась "
                f"с кодом {result.returncode}. Переход к следующей команде."
            )
        
        if result.returncode != 0:
            print(f"Текст ошибки записан в: {log_path}")

            with open(log_path, "a", encoding="utf-8") as log_file:
                log_file.write(f"[{datetime.now().isoformat(timespec='seconds')}] "
                                f"Команда: {' '.join(cmd)}\n")
                log_file.write(f"Код возврата: {result.returncode}\n")
                log_file.write("Stderr:\n")
                log_file.write(result.stderr or "(пусто)")
                log_file.write("\n\n")

    print("\nВсе команды выполнены.")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Запускает последовательность команд npge в указанной директории."
    )
    parser.add_argument(
        "path",
        help="Путь до папки, в которой нужно выполнить команды npge",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_commands(args.path)
