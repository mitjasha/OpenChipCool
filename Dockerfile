FROM python:3.12-slim
WORKDIR /work
COPY pyproject.toml README.md ./
COPY src ./src
COPY dashboard.py ./
COPY results ./results
RUN pip install --no-cache-dir '.[dashboard]'
EXPOSE 8501
CMD ["streamlit","run","dashboard.py","--server.address=0.0.0.0"]
