from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import io
import time
import numpy as np
import cv2
from PIL import Image
from tflite_runtime.interpreter import Interpreter

import os
from dotenv import load_dotenv
from flask import Flask, request, abort
from linebot import LineBotApi, WebhookHandler
from linebot.exceptions import InvalidSignatureError
from linebot.models import MessageEvent, TextMessage, TextSendMessage, QuickReply, QuickReplyButton, MessageAction

# 加載 .env 文件中的變數
load_dotenv()

# 從環境變數中讀取 LINE 的 Channel Access Token 和 Channel Secret
line_token = os.getenv('LINE_TOKEN')
line_secret = os.getenv('LINE_SECRET')

# 檢查是否設置了環境變數
if not line_token or not line_secret:
    print(f"LINE_TOKEN: {line_token}")  # 調試輸出
    print(f"LINE_SECRET: {line_secret}")  # 調試輸出
    raise ValueError("LINE_TOKEN 或 LINE_SECRET 未設置")

# 初始化 LineBotApi 和 WebhookHandler
line_bot_api = LineBotApi(line_token)
handler = WebhookHandler(line_secret)

# 創建 Flask 應用
app = Flask(__name__)

# 設置一個路由來處理 LINE Webhook 的回調請求
@app.route("/", methods=['POST'])
def callback():
    # 取得 X-Line-Signature 標頭
    signature = request.headers['X-Line-Signature']

    # 取得請求的原始內容
    body = request.get_data(as_text=True)
    app.logger.info(f"Request body: {body}")

    # 驗證簽名並處理請求
    try:
        handler.handle(body, signature)
    except InvalidSignatureError:
        abort(400)

    return 'OK'

# 自動發送訊息給用戶
def push_message_to_user(user_id):
    global total  # 使用全局變量
    try:
        # 發送訊息給用戶
        push_text = f"自動通知：目前的總結帳金額是: {total} dollars"
        line_bot_api.push_message(user_id, TextSendMessage(text=push_text))
        
        message = TextSendMessage(
            text='請點選想要的結帳方式~',
            quick_reply=QuickReply(
                items=[
                    QuickReplyButton(
                        action=MessageAction(label="現金", text="現金")
                    ),
                    QuickReplyButton(
                        action=MessageAction(label="Line_pay", text="Line_pay")
                    ),
                    QuickReplyButton(
                        action=MessageAction(label="信用卡", text="信用卡")
                    ),
                ]
            )
        )
        line_bot_api.push_message(user_id, message)
    except Exception as e:
        line_bot_api.push_message(user_id,
            TextSendMessage(text='Sorry~故障囉！'))
        print(f"Error: {e}")



def load_labels(path):
    with open(path, 'r') as f:
        return {i: line.strip() for i, line in enumerate(f.readlines())}

def set_input_tensor(interpreter, image):
    tensor_index = interpreter.get_input_details()[0]['index']
    input_tensor = interpreter.tensor(tensor_index)()[0]
    input_tensor[:, :] = image

def classify_image(interpreter, image, top_k=1):
    """Returns a sorted array of classification results."""
    set_input_tensor(interpreter, image)
    interpreter.invoke()
    output_details = interpreter.get_output_details()[0]
    output = np.squeeze(interpreter.get_tensor(output_details['index']))

    # If the model is quantized (uint8 data), then dequantize the results
    if output_details['dtype'] == np.uint8:
        scale, zero_point = output_details['quantization']
        output = scale * (output - zero_point)

    ordered = np.argpartition(-output, top_k)
    return [(i, output[i]) for i in ordered[:top_k]]

def main():
    global total  # 使用全局變量
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--model', help='File path of .tflite file.', required=True)
    parser.add_argument(
        '--labels', help='File path of labels file.', required=True)
    args = parser.parse_args()

    labels = load_labels(args.labels)

    #interpreter = tf.lite.Interpreter(args.model)
    interpreter = Interpreter(args.model)

    interpreter.allocate_tensors()
    _, height, width, _ = interpreter.get_input_details()[0]['shape']

    #with picamera.PiCamera(resolution=(640, 480), framerate=30) as camera:
        #camera.start_preview()
    cap = cv2.VideoCapture(0)
    #cap = cv2.VideoCapture(0,cv2.CAP_GSTREAMER)
    #cap = cv2.VideoCapture(0,cv2.CAP_V4L)
    #擷取畫面 寬度 設定為640
    cap.set(cv2.CAP_PROP_FRAME_WIDTH,640)
    #擷取畫面 高度 設定為480
    cap.set(cv2.CAP_PROP_FRAME_HEIGHT, 480)

    key_detect = 0
    times = 1
    last_label_id = None
    label_start_time = None
    total = 0
    price = 0
    price_added = False  # 新增標誌變量

    # 定義標籤和價格的對應關係
    price_dict = {
        "None": 0,
        "pen": 10,
        "drink": 30,
        "cookie": 8,
        "shampoo": 50,
        "tissue": 20
    }

    while (key_detect == 0):
        ret, image_src = cap.read()
        print(image_src)
        frame_width = image_src.shape[1]
        frame_height = image_src.shape[0]

        cut_d = int((frame_width - frame_height) / 2)
        crop_img = image_src[0:frame_height, cut_d:(cut_d + frame_height)]

        image = cv2.resize(crop_img, (224, 224), interpolation=cv2.INTER_AREA)

        start_time = time.time()
        if (times == 1):
            results = classify_image(interpreter, image)
            elapsed_ms = (time.time() - start_time) * 1000
            label_id, prob = results[0]
            if label_id == last_label_id:
                if label_start_time is None:
                    label_start_time = time.time()
                elif time.time() - label_start_time >= 2.5:
                    print("Product" + labels[label_id], prob)
                    if not price_added:  # 檢查是否已經累加過
                        price = price_dict.get(labels[label_id], 0)  # 使用字典獲取價格，默認為0
                        if price != 0 :
                            sub_price = price 
                        total = total + price
                        price_added = True  # 設置標誌變量為 True
            else:
                last_label_id = label_id
                label_start_time = time.time()
                price_added = False  # 重置標誌變量

        # 在每次循環中都顯示文本
        cv2.putText(crop_img, "Product : " + labels[label_id] + " " + str(price) + " dollars", (5, 30), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 1, cv2.LINE_AA)
        cv2.putText(crop_img, "Total : " + str(total) + " dollars", (5, 60), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 1, cv2.LINE_AA)
        
        times = times + 1
        if (times > 1):
            times = 1
        
        key_pressed = cv2.waitKey(1) & 0xFF  # 只檢測一次按鍵

        if key_pressed == ord('-'):
            total -= sub_price

        elif key_pressed == ord('r'):  # 重置
            total = 0
            label_id = None
            price_added = False  # 重置標誌變量

        elif key_pressed == ord('q'):
            key_detect = 1

        elif key_pressed == ord('c'):  # 按 c 鍵推送訊息
            user_id = 'U3eacbf287900d372524b2c95b544af6a'  # 替換為實際的用戶 ID
            push_message_to_user(user_id)  # 傳遞 total 參數

        
        cv2.imshow('Detecting....', crop_img)

    cap.release()
    cv2.destroyAllWindows()

if __name__ == '__main__':
    main()
    app.run(host='0.0.0.0', port=5000)
