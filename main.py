from fastapi import FastAPI, UploadFile, File
import shutil
import subprocess
import os
import traceback

app = FastAPI()

@app.post("/nest/")
async def nest_file(file: UploadFile = File(...)):
    try:
        # Создание директорий для файлов
        os.makedirs("uploads", exist_ok=True)
        os.makedirs("outputs", exist_ok=True)

        # Пути input и output
        input_path = f"uploads/{file.filename}"
        output_path = f"outputs/nested_{file.filename}"

        # Сохраняем загруженный файл
        with open(input_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)

        # Запускаем Nesting CLI
        nest_binary = "nest.exe" if os.name == "nt" else "./nest"
        result = subprocess.run(
            [nest_binary, input_path, output_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Проверка на ошибки
        if result.returncode != 0:
            return {
                "error": result.stderr,
                "stdout": result.stdout,
                "returncode": result.returncode
            }

        return {
            "output_file": output_path,
            "stdout": result.stdout
        }

    except Exception as e:
        return {
            "exception": str(e),
            "traceback": traceback.format_exc()
        }
