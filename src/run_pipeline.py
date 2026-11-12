import argparse
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import seq_auto

SCRIPT_DIR = Path(__file__).resolve().parent
SEQ_AUTO_SCRIPT = SCRIPT_DIR / "seq_auto.py"
AUTO_NPGE_SCRIPT = SCRIPT_DIR / "auto_npge.py"
ALAN_SCRIPT = SCRIPT_DIR / "alan.py"

ALAN_OUTDIR_NAME = "alan_result"

LOG_NAME = "pipeline_errors.log"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Последовательно запускает seq_auto.py (сбор последовательностей по "
            "таблице сборок из NCBI), auto_npge.py (построение пангенома NPGE) "
            "и alan.py (анализ результатов) в указанной рабочей директории. "
            "Флаги -s/-n/-a позволяют запустить только один шаг. Каждый шаг выполняется отдельным "
            f"процессом, ошибки записываются в {LOG_NAME} внутри рабочей директории."
        )
    )
    parser.add_argument(
            "directory",
            help=(
                "Рабочая директория NPGE: в неё сохраняется файл с последовательностями "
                "и в ней же выполняются команды npge"
            ),
    )
    parser.add_argument(
        "-t", "--table",
        default=None,
        help=(
            "Путь к входной таблице сборок (.tsv или .csv) для seq_auto "
            "(обязателен, если запускается seq_auto)"
        ),
    )

    # Выбор шагов (если не указан ни один — выполняются все)
    parser.add_argument(
        "-s", "--only-seq-auto",
        action="store_true",
        help="Запустить только seq_auto",
    )
    parser.add_argument(
        "-n", "--only-npge",
        action="store_true",
        help="Запустить только auto_npge",
    )
    parser.add_argument(
        "-a", "--only-alan",
        action="store_true",
        help="Запустить только alan",
    )

    # Параметры auto_npge
    parser.add_argument("-i", "--identity", help = "Минимальное сходство для блока, доля")
    parser.add_argument("-l", "--length", help = "Минимальная длина блока, пары нуклеотидов")

    # Параметры alan
    parser.add_argument(
        "-c", "--csv",
        action="store_true",
        help="alan: сохранять таблицы в csv вместо Excel (флаг -c в alan.py)",
    )

    # Параметры seq_auto
    parser.add_argument(
        "-e", "--email",
        default="your_email@example.com",
        help="Email для Entrez.email (обязателен для запросов к NCBI)",
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="NCBI API-ключ (опционально, увеличивает лимит запросов)",
    )
    parser.add_argument(
        "--assembly-col",
        default=seq_auto.ASSEMBLY_COL,
        help=f"Название столбца с accession сборки (по умолчанию: '{seq_auto.ASSEMBLY_COL}')",
    )
    parser.add_argument(
        "--linkname",
        default=seq_auto.LINKNAME,
        help=f"Имя elink-связи между assembly и nucleotide (по умолчанию: '{seq_auto.LINKNAME}')",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=seq_auto.MAX_RETRIES,
        help=f"Максимальное число повторных попыток Entrez-запроса (по умолчанию: {seq_auto.MAX_RETRIES})",
    )
    parser.add_argument(
        "--base-delay",
        type=float,
        default=seq_auto.BASE_DELAY,
        help=f"Базовая задержка перед повтором запроса, секунды (по умолчанию: {seq_auto.BASE_DELAY})",
    )
    parser.add_argument(
        "--request-delay",
        type=float,
        default=seq_auto.REQUEST_DELAY,
        help=f"Задержка между запросами к NCBI, секунды (по умолчанию: {seq_auto.REQUEST_DELAY})",
    )

    args = parser.parse_args()

    if "seq_auto" in selected_steps(args) and not args.table:
        parser.error("для запуска seq_auto необходимо указать таблицу сборок: -t/--table")

    return args


def selected_steps(args) -> list:
    """
    Возвращает имена шагов, которые нужно выполнить, в порядке запуска.
    Если ни один из флагов -s/-n/-a не указан — выполняются все шаги.
    """
    run_all = not (args.only_seq_auto or args.only_npge or args.only_alan)
    steps = []
    if run_all or args.only_seq_auto:
        steps.append("seq_auto")
    if run_all or args.only_npge:
        steps.append("auto_npge")
    if run_all or args.only_alan:
        steps.append("alan")
    return steps


def log_failure(log_path: Path, cmd: list, returncode, stderr: str) -> None:
    """Дописывает информацию об упавшей команде в лог-файл."""
    with open(log_path, "a", encoding="utf-8") as log_file:
        log_file.write(f"[{datetime.now().isoformat(timespec='seconds')}] "
                       f"Команда: {' '.join(map(str, cmd))}\n")
        log_file.write(f"Код возврата: {returncode}\n")
        log_file.write("Stderr:\n")
        log_file.write(stderr.strip() if stderr and stderr.strip() else "(пусто)")
        log_file.write("\n\n")


def run_step(cmd: list, log_path: Path, cwd: Path) -> bool:
    """
    Запускает команду как отдельный процесс. stdout идёт прямо в консоль
    (виден прогресс), stderr перехватывается. При ненулевом коде возврата
    или невозможности запустить процесс ошибка пишется в лог.
    Возвращает True, если шаг завершился успешно.
    """
    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            env=env,
            text=True,
            encoding="utf-8",
            errors="replace",
            stderr=subprocess.PIPE,
        )
    except OSError as e:
        print(f"Не удалось запустить команду: {e}")
        log_failure(log_path, cmd, "не запущена", str(e))
        print(f"Текст ошибки записан в: {log_path}")
        return False

    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")

    if result.returncode != 0:
        print(f"\nКоманда завершилась с кодом {result.returncode}. "
              f"Текст ошибки записан в: {log_path}")
        log_failure(log_path, cmd, result.returncode, result.stderr)
        return False

    return True


def main():
    args = parse_args()

    directory = Path(args.directory).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    log_path = directory / LOG_NAME

    steps = selected_steps(args)
    commands = {}

    if "seq_auto" in steps:
        seq_cmd = [
            sys.executable, str(SEQ_AUTO_SCRIPT),
            str(Path(args.table).expanduser().resolve()),
            str(directory),
            "-e", args.email,
            "--assembly-col", args.assembly_col,
            "--linkname", args.linkname,
            "--max-retries", str(args.max_retries),
            "--base-delay", str(args.base_delay),
            "--request-delay", str(args.request_delay),
        ]
        if args.api_key:
            seq_cmd += ["--api-key", args.api_key]
        commands["seq_auto"] = seq_cmd

    commands["auto_npge"] = [sys.executable, str(AUTO_NPGE_SCRIPT), str(directory)]
    if args.identity is not None:
        commands["auto_npge"] += ["-i", str(args.identity)]
    if args.length is not None:
        commands["auto_npge"] += ["-l", str(args.length)]

    commands["alan"] = [
        sys.executable, str(ALAN_SCRIPT),
        str(directory),
        str(directory / ALAN_OUTDIR_NAME),
    ]
    if args.csv:
        commands["alan"].append("-c")

    for number, name in enumerate(steps, start=1):
        print(f"\n=== Шаг {number}/{len(steps)}: {name} ===")
        if not run_step(commands[name], log_path, directory):
            print(f"{name} завершился с ошибкой — дальнейшие шаги не выполняются.")
            sys.exit(1)

    print("\nПайплайн завершён.")


if __name__ == "__main__":
    main()
